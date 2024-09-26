"""
QTL calling
"""

import logging
import dataclasses
from dataclasses import dataclass
from typing import TypedDict

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.special import betainc

from . import vcf
from . import stats

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def split_genes(expression_df: pd.DataFrame, n_bins: int
        ) -> list[tuple[int, int]]:
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

    return _filter_genotype_sites(genotype, keep_sites)

def _filter_genotype_sites(genotype, keep_sites):
    return _Genotype(
            genotype.samples,
            genotype.meta[keep_sites],
            genotype.valid[keep_sites, :],
            genotype.dosages[keep_sites, :],
            )

@dataclass(frozen=True)
class _Genotype:
    samples: list[str]
    meta: pd.DataFrame
    valid: npt.NDArray[np.bool]
    dosages: npt.NDArray[np.uint8]

class QtlConfigDict(TypedDict, total=False):
    """
    Qtl calling options, see :func:`call_qtls`
    """
    window_size: int
    maf: float
    impute_genotype: bool
    num_null_genes: int

DEFAULT_QTL_CONFIG: QtlConfigDict = {
        'window_size': 10**6,
        'maf': 0.01,
        'impute_genotype': True,
        'num_null_genes': 10**4,
        }
"""
Default options for :func:`call_qtls`
"""

def call_qtls(expression_df: tuple[pd.DataFrame, int], gene_window_indexes:
        tuple[int, int], vcf_path: str, gt_covariates_df: pd.DataFrame,
        gx_covariates_df: pd.DataFrame | None,
        vcf_index: vcf.VCFSimpleIndex | None = None,
        qtl_config: QtlConfigDict | None = None):
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
        gt_covariates_df: covariate data frame used to residualize genotype.
            First column is interpreted as covariate identifier, all other
            columns as values for the sample of matching column name.
        gx_covariates_df: different set of covariates to use specifically for
            gene expression. Can be None for no covariates.
        vcf_index: simple index generated by :func:`gtex_qtl.vcf.make_simple_index`.
            If not given, will be generated on the fly.
        qtl_config: dictionnary of QTL calling options

            - **window_size**: number of bases upstream and downstream of the `start`
            - **maf**: minimum maf filter to apply
            - **impute_genotype**: whether to replace missing genotype dosages
              by the sample mean.
            - **num_null_genes**: number of genes to sample to create the null
              set for empirical null distribution estimation

            Defaults taken from :data:`DEFAULT_QTL_CONFIG`

    Returns:
        data frame of pair associations
    """
    merged_qtl_config = DEFAULT_QTL_CONFIG | (qtl_config or {})

    expression = _pack_expression(expression_df)

    chrom, genotype_window = _make_genotype_window(expression['meta_x'],
            gene_window_indexes, merged_qtl_config['window_size'])

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
    genotype = _filter_genotype_maf(genotype, merged_qtl_config['maf'])
    genotype = _filter_genotype_samples(genotype, expression['samples'])

    if merged_qtl_config['impute_genotype']:
        logger.info('Imputing missing genotypes')
        genotype = _impute_genotypes(genotype)

    gt_covariates = _pack_covariates(gt_covariates_df, expression['samples'])
    logger.info('Loaded %d genotype covariates', len(gt_covariates['meta_c']))

    if gx_covariates_df is None:
        gx_covariates = None
    else:
        gx_covariates = _pack_covariates(gx_covariates_df, expression['samples'])
        logger.info('Loaded %d gene expression covariates', len(gx_covariates['meta_c']))


    genotype = _regress_genotype(genotype, gt_covariates)

    logger.info('Computing associations')
    pairs, summaries_perm, summaries_ic = zip(*(
        _call_gene(genotype, expression_item, expression, gx_covariates, merged_qtl_config)
        for expression_item in _iter_expression(expression, gene_window_indexes)
        ))

    return (
        # At this point the pairs are still indexed by the original index of the
        # sites, which is both sparse and contain duplicates, so reset it.
        pd.concat(pairs).reset_index(drop=True),
        pd.DataFrame(summaries_perm), pd.DataFrame(summaries_ic)
        )

class _CovariatesDict(TypedDict):
    """
    Split of covariate table into metadata table and values array
    """
    meta_c: pd.DataFrame
    values_cs: npt.NDArray[np.float32]

def _pack_covariates(covariates_df: pd.DataFrame, expression_samples: list[str]
        ) -> _CovariatesDict:
    """
    Split covariate df into a metadata frame and a compact array of values,
    after subsetting the right columns
    """
    return {
            'meta_c': covariates_df[[covariates_df.columns[0]]],
            'values_cs': np.array(covariates_df[expression_samples])
            }

class _ExpressionDict(TypedDict):
    """
    Expression dataframe with expression values extracted as an array
    """
    meta_x: pd.DataFrame
    values_xs: npt.NDArray[np.float32]
    samples: list[str]

def _pack_expression(expression: tuple[pd.DataFrame, int]) -> _ExpressionDict:
    """
    Split expression data frame
    """
    expression_df, n_meta = expression

    meta_columns = expression_df.columns[:n_meta]
    data_columns = expression_df.columns[n_meta:]

    return {
            'meta_x': expression_df[meta_columns],
            'values_xs': np.array(expression_df[data_columns]),
            'samples': list(data_columns)
            }

def _regress_genotype(genotype, covariates):
    return dataclasses.replace(genotype,
            dosages=stats.regress_missing(
                genotype.dosages,
                genotype.valid,
                covariates['values_cs']
                )
            )

def _regress_expression(expression, covariates):
    if covariates is None:
        return expression
    return dict(expression,
            values_xs=stats.regress(
                expression['values_xs'],
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
            expression['meta_x'].iloc[start:stop].iterrows(),
            expression['values_xs'][start:stop]
            )
        )

def _call_pairs(genotype, expression):
    gt_gs = genotype.dosages
    valid_gs = genotype.valid

    gx_xs = expression['values_xs']

    gtgx_xg = gx_xs @ gt_gs.T
    gt2_g = np.sum(gt_gs**2, -1)

    slope_xg = gtgx_xg / gt2_g

    res2_xg = (
            (gx_xs**2) @ valid_gs.T
            - 2. * gtgx_xg * slope_xg
            + slope_xg**2 * gt2_g
            )

    # Scaled t^2 statistic (meaning up to the dofs factor)
    st2_xg = slope_xg**2 * gt2_g / res2_xg

    return slope_xg, st2_xg

def _dev_and_pval(slope_xg, st2_xg, n_valid_g, dofs):
    """
    Estimate standard deviation and p-value for a given number of
    missing degrees of freedom
    """

    # Residual degrees of freedom: number of valid observations, minus the
    # parameters eliminated by pre-regression -- including the mean --, minus
    # the genotype weight parameter itself
    rdofs_g = n_valid_g - dofs - 1

    slope_var_xg = slope_xg**2 / (rdofs_g * st2_xg)

    slope_pval_xg = betainc(.5 * rdofs_g, .5, 1. / (1. + st2_xg))

    return slope_pval_xg, np.sqrt(slope_var_xg)

def _make_null_genes(expression_item, n_perms):
    """
    Generate expression with no association to genotype

    Return an expression dictionnary looking like the expression of several
    genes filled with the null set of genes (or decoys).
    """

    ## Seed rng with pos, so that permutation is different for each gene but
    ## replicable
    rng = np.random.default_rng(expression_item['meta']['start'])

    values_s = expression_item['values_s']

    n_samples = len(values_s)
    # Permute by drawing random reals and sorting them to get random orders
    perm_values_xs = values_s[
            np.argsort(
                rng.random((n_perms, n_samples)),
                -1)
            ]

    return dict(expression_item,
            values_xs=perm_values_xs)

def _filter_genotype_proximity(genotype, expression_item, qtl_config):
    """
    Filter genotype sites within requested window of gene feature start
    """
    rel_pos = genotype.meta['POS'] - expression_item['meta']['start']

    genotype = _filter_genotype_sites(genotype,
            np.abs(rel_pos) <= qtl_config['window_size']
            )
    return  dataclasses.replace(genotype,
            # Not actually a distance
            meta=genotype.meta.assign(tss_distance=rel_pos)
            )

def _call_gene(genotype, expression_item, expression, covariates,
        qtl_config: QtlConfigDict):
    """
    Compute all pairwise associations, plus the gene-level statistics for one
    gene

    Args:
        genotype: all sites to consider, needs not be filtered by position yet
    """

    # Filter genotype by distance
    genotype = _filter_genotype_proximity(genotype, expression_item, qtl_config)

    pairs = _call_true_pairs(genotype, expression_item, covariates)

    best_pair = pairs.iloc[pairs.pval_nominal.argmin()]

    summary_perm = _call_null_perm(
            genotype, expression_item, covariates, qtl_config,
            best_pair
            )

    summary_interchrom = _call_null_interchrom(
            genotype, expression_item, expression,
            covariates, qtl_config, best_pair
            )

    return pairs, summary_perm, summary_interchrom

def _call_true_pairs(genotype, expression_item, covariates):
    """
    Compute association on the real gene expression and generate dataframe of
    associations
    """

    # Make the expression look like an array of genes with one gene
    true_expression = dict(expression_item,
            values_xs=expression_item['values_s']
            )
    del true_expression['values_s']

    # Residualize true gene
    true_expression = _regress_expression(true_expression, covariates)

    # Compute associations
    slope_g, scaled_t2_g = _call_pairs(genotype, true_expression)

    # Compute deviation and tests ("nominal" p-values)

    ## Dofs is number of covariates plus 1 since there's an intercept
    if covariates is None:
        dofs = 0
    else:
        dofs = len(covariates['meta_c']) + 1

    slope_pval_g, slope_std_g = _dev_and_pval(
            slope_g, scaled_t2_g, np.sum(genotype.valid, -1), dofs
            )

    # Gather associations in a dataframe
    return (
        genotype.meta
        .assign(**expression_item['meta'])
        .assign(
            pval_nominal=slope_pval_g,
            slope=slope_g,
            slope_se=slope_std_g,
            slope_st2=scaled_t2_g,
            )
        )

def _call_null_perm(genotype, expression_item, covariates, qtl_config, best_pair):
    """
    Compute null distribution statistics from permutations
    """

    # Generate decoy genes
    null_expression = _make_null_genes(expression_item,
            qtl_config['num_null_genes'])

    # Residualize decoy genes
    null_expression = _regress_expression(null_expression, covariates)

    return _call_null_any(genotype, expression_item, null_expression, best_pair)

def _call_null_any(genotype, expression_item, null_expression, best_pair):
    """
    Compute null distribution statistics from a pre-defined set of null genes
    """
    # Compute associations (only the test statistic matters)
    _, null_scaled_t2_xg = _call_pairs(genotype, null_expression)

    # Estimate dofs from decoy genes
    ## Keep best of each decoy
    _best_null_st2s = np.max(null_scaled_t2_xg, -1)

    ## Estimate both meaningful parameters jointly
    jmle_params = stats.fit_max_scaled_t(_best_null_st2s)
    # Same with correct formula
    jmle_params_rev2 = stats.fit_max_scaled_t(_best_null_st2s, revision=2)

    ## Estimate parameters in turn with fqtl procedure
    fqtl_params = stats.fqtl_fit_max_scaled_t(_best_null_st2s)

    return dict(expression_item['meta'],
            num_var=len(genotype.meta),
            # FQTL reproduction
            **fqtl_params,
            **stats.fqtl_pval_max_scaled_t(best_pair.slope_st2, fqtl_params),
            # Joint MLE estimation
            **jmle_params,
            **{ k + '_jmle': v for k, v in
                stats.pval_max_scaled_t(best_pair.slope_st2, jmle_params).items()},
            # Revision
            **{ k + '_rev2': v for k, v in jmle_params_rev2.items()},
            **{ k + '_jmle_rev2': v for k, v in
                stats.pval_max_scaled_t(best_pair.slope_st2, jmle_params_rev2).items()},
            # Raw data and permutation based
            **{k: best_pair[k]
                for k in ('pval_nominal', 'slope', 'slope_se')
                },
            pval_perm=_empirical_pval(best_pair.slope_st2, _best_null_st2s)
            )


def _call_null_interchrom(genotype, expression_item, expression, covariates,
        qtl_config, best_pair):
    """
    Compute null distribution statistics from between-chromosome associations
    """

    null_expression = _choose_interchrom_genes(expression_item, expression,
            qtl_config['num_null_genes'])

    # Residualize decoy genes
    # In the future, this would be done ahead of this function
    null_expression = _regress_expression(null_expression, covariates)

    return _call_null_any(genotype, expression_item, null_expression, best_pair)

def _choose_interchrom_genes(expression_item, expression, num_genes):
    """
    Sample genes from a different chromosome
    """
    # BED convention
    chr_col = expression['meta_x'].columns[0]

    target_chr = expression_item['meta'][chr_col]

    igenes_indices, = np.nonzero(np.array(
        expression['meta_x'][chr_col] != target_chr
        ))

    logger.info('Found %d inter-chrom genes to choose from out of %d',
            len(igenes_indices), len(expression['meta_x']))

    # Same rng seeding strategy as in _make_null_genes
    rng = np.random.default_rng(expression_item['meta']['start'])
    null_indices = rng.choice(igenes_indices, num_genes, replace=False)

    return {'values_xs': expression['values_xs'][null_indices]}

def _empirical_pval(stat, null_stats):
    """
    Empirical p-value given samples of the test statistic

    P-value is computed for the upper tail of the statistic.
    """

    return (
            (np.sum(null_stats >= stat) + 1)
            / (len(null_stats) + 1)
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
