"""
QTL calling
"""

import logging
from dataclasses import dataclass
from typing import List

import numpy as np

from . import vcf

logger = logging.getLogger(__name__)

def split_genes(expression_df, n_bins):
    """
    Allocates genes in bins such that bins are within the same chromosomes

    Args:
        expression_df: Pandas data frame of genes with the chromosome id as the
            first column (this is in line with the BED file format, which has
            columns standardized by position but not by name) n_bins: number of
            bins to make. If less than the number of chromosomes
            in `expression_df`, it will be increased to that number.

    Returns: list of pairs of indexes, corresponding to the first and last gene
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
    mafs = genotype.dosages.sum(-1) / genotype.valid.sum(-1)
    mafs = np.minimum(mafs, 1. - mafs)

    keep_sites = mafs >= maf

    logger.info('Excluding %s sites with MAF less than %s '
            '(%s [%.2f%%] remaining)',
            len(keep_sites) - keep_sites.sum(), maf, keep_sites.sum(),
            100 * keep_sites.sum() / len(keep_sites)
            )

    filtered = _Genotype(
            genotype.samples,
            [rec for rec, keep_site in zip(genotype.meta, keep_sites)],
            genotype.valid[keep_sites, :],
            genotype.dosages[keep_sites, :],
            )

    return filtered

@dataclass
class _Genotype:
    samples: List[str]
    meta: List[tuple]
    valid: np.array
    dosages: np.array

def call_qtls(expression_df, gene_window_indexes, vcf_path,
        vcf_index=None, window_size=10**6, maf=0.01):
    """
    Draft QTL calling.

    Args:
        expression_df: gene expression data frame. First column must be the
            chromosome id, second the start gene of the gene feature (`start`
            thereafter) ; all other columns are phenotypes.
        gene_window_indexes: tuple (first, last) of numeric indexes, inclusive,
            delimiting the rows of `expression_df` to process. Both rows, and
            all the ones in between must be located in the same chromosome and
            sorted by the `start` field.
        window_size: number of bases upstream and downstream of the `start`
        vcf_path: path to the BGZIP-compressed vcf containing the genotypes in
            "GT" format.
        vcf_index: simple index generated by :func:`gtex_qtl.vcf.make_simple_index`.
            If not given, will be generated on the fly.
    """
    chrom, genotype_window = _make_genotype_window(expression_df, gene_window_indexes, window_size)

    _, _, *expression_samples = expression_df.columns
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
    genotype = _filter_genotype_maf(genotype, maf)
    genotype = _filter_genotype_samples(genotype, expression_samples)
