"""
QTL calling
"""

import logging
import dataclasses
from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd
from scipy.special import erfc

from . import vcf
from . import stats

logger = logging.getLogger(__name__)

def split_genes(expression_df, n_bins):
    """
    Allocates genes in bins such that bins are within the same chromosomes

    Args:
        expression_df: Pandas data frame of genes with the chromosome id as the
            first column (this is in line with the BED file format, which has
            columns standardized by position but not by name)
        n_bins: number of
            bins to make. If less than the number of chromosomes
            in `expression_df`, it will be increased to that number.

    Returns:
        list of pairs of indexes, corresponding to the first and last gene
        of each bin, inclusive.
    """

    chrom_col = expression_df.columns[0]

    grouped_edf = expression_df.reset_index().groupby(chrom_col)
    genes_per_chrom = grouped_edf.size()

    n_bins = max(n_bins, len(genes_per_chrom))

    # Allocate one bin to each chrom
    bins = genes_per_chrom // genes_per_chrom

    # Allocate the rest priorizing minimizing max bin size
    for _ in range(bins.sum(), n_bins):
        largest_bin = (genes_per_chrom / bins).argmax()
        bins.iloc[largest_bin] += 1

    assert bins.sum() == n_bins

    # Finally, generate the required number of intervals per chrom
    chrom_bounds = (
            grouped_edf['index']
            .agg([np.min])
            .assign(bins=bins, genes=genes_per_chrom)
            )


    # Mind the rounding: int(x + .5) will consistently give the integer x is
    # close to and int(x - .5) the integer below when x is close to an integer
    # boundary
    pairs = [(
            int(start + i * genes / c_bins + .5),
            int(start + (i+1) * genes / c_bins - .5)
            )
        for _, (start, c_bins, genes) in chrom_bounds.iterrows()
        for i in range(c_bins)
        ]

    assert sum(last - first + 1 for first, last in pairs) == len(expression_df)
    assert all(expression_df[chrom_col].iloc[first] == expression_df[chrom_col].iloc[last]
            for first, last in pairs)

    return pairs

def _make_genotype_window(expression_df, gene_window_indexes, window_size):
    chrom_col, start_col = expression_df.columns[:2]

    gene_meta = expression_df[[chrom_col, start_col]]
    (_, first_meta), (_, last_meta) = (
        gene_meta
        .iloc[list(gene_window_indexes)]
        .iterrows()
        )

    chrom, chrom_last = first_meta[chrom_col], last_meta[chrom_col]

    if chrom != chrom_last:
        raise ValueError('First and last gene must be in the same chromosome '
            f'(got "{chrom}" and "{chrom_last}").')

    return chrom, (
            int(first_meta[start_col]) - window_size,
            int(last_meta[start_col]) + window_size
            )

def _filter_genotype_samples(genotype, expression_samples):
    esample_set = set(expression_samples)
    keep_samples = np.array([
        sample in esample_set
        for sample in genotype.samples])
    logger.info('Excluding %s samples with no matching phenotype from genotype '
            '(%s remaining)',
            len(genotype.samples) - keep_samples.sum(), keep_samples.sum()
            )

    filtered = _Genotype(
            [sample for sample in genotype.samples if sample in esample_set],
            genotype.meta,
            genotype.valid[:, keep_samples],
            genotype.dosages[:, keep_samples]
            )

    return filtered

def _filter_genotype_maf(genotype, maf):
    mafs = genotype.dosages.sum(-1) / genotype.valid.sum(-1) / 2.
    mafs = np.minimum(mafs, 1. - mafs)

    keep_sites = mafs >= maf

    logger.info('Excluding %s sites with MAF less than %s '
            '(%s [%.2f%%] remaining)',
            len(keep_sites) - keep_sites.sum(), maf, keep_sites.sum(),
            100 * keep_sites.sum() / len(keep_sites)
            )

    filtered = _Genotype(
            genotype.samples,
            genotype.meta[keep_sites],
            genotype.valid[keep_sites, :],
            genotype.dosages[keep_sites, :],
            )

    return filtered

@dataclass(frozen=True)
class _Genotype:
    samples: List[str]
    meta: pd.DataFrame
    valid: np.array
    dosages: np.array

DEFAULT_QTL_CONFIG = {
        'window_size': 10**6,
        'maf': 0.01,
        'impute_genotype': False
        }
"""
Default options for :func:`call_qtls`
"""

def call_qtls(expression_df, gene_window_indexes, vcf_path,
        covariates_df, vcf_index=None, qtl_config=None):
    """
    Draft QTL calling.

    Args:
        expression_df: tuple with gene expression data frame and number of
            initial metadata columns. First column must be the
            chromosome id, second the start gene of the gene feature (`start`
            thereafter) ; all other columns after the number of given metadata
            columns are phenotypes.
        gene_window_indexes: tuple (first, last) of numeric indexes, 0-based,
            inclusive, delimiting the rows of `expression_df` to process. Both
            rows, and all the ones in between must be located in the same
            chromosome and sorted by the `start` field.
        vcf_path: path to the BGZIP-compressed vcf containing the genotypes in
            "GT" format.
        covariates_df: covariate data frame. First column is interpreted as
            covariate identifier, all other columns as values for the sample of
            mathching column name.
        vcf_index: simple index generated by :func:`gtex_qtl.vcf.make_simple_index`.
            If not given, will be generated on the fly.
        qtl_config: dictionnary of QTL calling options

            - **window_size**: number of bases upstream and downstream of the `start`
            - **maf**: minimum maf filter to apply
            - **impute_genotype**: whether to replace missing genotype dosages
              by the sample mean.

            Defaults taken from :data:`DEFAULT_QTL_CONFIG`

    Returns:
        data frame of pair associations
    """
    qtl_config = dict(DEFAULT_QTL_CONFIG, **(qtl_config or {}))

    expression = _pack_expression(expression_df)

    chrom, genotype_window = _make_genotype_window(expression['meta_g'],
            gene_window_indexes, qtl_config['window_size'])

    _, genotype_samples = vcf.read_header(vcf_path)
    if vcf_index is None:
        logger.info('Indexing VCF')
        vcf_index = vcf.make_simple_index(vcf_path)

    logger.info('Reading genotypes from VCF')
    genotype = _Genotype(genotype_samples, *vcf.parse_region_indexed(
                vcf_path, vcf_index,
                chrom, *genotype_window)
                )

    logger.info("Loaded %s sites from %s samples", genotype.dosages.shape[0],
            genotype.dosages.shape[1])
    logger.info("%.3f%% missing genotypes",
            100 * (1. - np.sum(genotype.valid) / genotype.valid.size))

    # Mind the order !
    genotype = _filter_genotype_maf(genotype, qtl_config['maf'])
    genotype = _filter_genotype_samples(genotype, expression['samples'])

    if qtl_config['impute_genotype']:
        logger.info('Imputing missing genotypes')
        genotype = _impute_genotypes(genotype)

    covariates = _pack_covariates(covariates_df, expression['samples'])
    genotype = _regress_genotype(genotype, covariates)
    expression = _regress_expression(expression, covariates)

    return pd.concat([
        # Dofs is number of covariates plus 1 since there's an intercept
        _compute_pairs(genotype, expression_item, len(covariates['meta_c']) + 1)
        for expression_item in _iter_expression(expression, gene_window_indexes)
        ])

def _pack_covariates(covariates_df, expression_samples):
    """
    Split covariate df into a metadata frame and a compact array of values,
    after subsetting the right columns
    """
    return dict(
        meta_c = covariates_df[[covariates_df.columns[0]]],
        values_cs = np.array(covariates_df[expression_samples])
        )

def _pack_expression(expression):
    """
    Split expression data frame
    """
    expression_df, n_meta = expression

    meta_columns = expression_df.columns[:n_meta]
    data_columns = expression_df.columns[n_meta:]

    return dict(
        meta_g = expression_df[meta_columns],
        values_gs = np.array(expression_df[data_columns]),
        samples = data_columns
        )

def _regress_genotype(genotype, covariates):
    return dataclasses.replace(genotype,
            dosages=stats.regress_missing(
                genotype.dosages,
                genotype.valid,
                covariates['values_cs']
                )
            )

def _regress_expression(expression, covariates):
    return dict(expression,
            values_gs=stats.regress(
                expression['values_gs'],
                covariates['values_cs']
                )
            )

def _iter_expression(expression, gene_window_indexes):
    start, stop = gene_window_indexes
    stop += 1
    return ({
        'meta': meta,
        'values_s': values,
        'samples': expression['samples']
        }
        for (_, meta), values in
        zip(
            expression['meta_g'].iloc[start:stop].iterrows(),
            expression['values_gs'][start:stop]
            )
        )

def _compute_pairs(genotype, expression_item, dofs):
    # FIXME: filter by position

    gt_gs = genotype.dosages
    valid_gs = genotype.valid
    n_valid_g = np.sum(valid_gs, -1)

    gx_s = expression_item['values_s']

    gtgx_g = gt_gs @ gx_s
    gt2_g = np.sum(gt_gs**2, -1)

    slope_g = gtgx_g / gt2_g

    res2_g = np.sum((gx_s - slope_g[:, None] * gt_gs)**2 * valid_gs, -1)

    slope_std_g = np.sqrt(res2_g  / (gt2_g * (n_valid_g - dofs - 1)))

    slope_t_g = slope_g / slope_std_g

    slope_normal_pval_g = erfc(np.abs(slope_t_g) / np.sqrt(2))

    return (
        genotype.meta
        .assign(**expression_item['meta'])
        .assign(
            slope=slope_g,
            slope_std=slope_std_g,
            slope_normal_pval=slope_normal_pval_g
            )
        )

def _impute_genotypes(genotype):
    """
    Replace missing genotype dosages by the sample mean
    """
    mean_g = np.sum(genotype.dosages, -1) / np.sum(genotype.valid, -1)

    return dataclasses.replace(genotype,
            dosages = np.where(genotype.valid, genotype.dosages, mean_g[:, None]),
            valid = genotype.valid |~ genotype.valid
            )
