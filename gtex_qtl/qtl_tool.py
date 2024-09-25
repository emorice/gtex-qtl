"""
Wrapper around the qtl caller for pipelining
"""

from typing import TypedDict

import pandas as pd
import galp
from galp import step

from . import qtl

class QtlToolConfigDict(TypedDict, total=False):
    """
    Options for qtl caller, see `call_qtl`
    """
    n_bins: int
    bed_fixed_fields: int
    qtl_core_config: qtl.QtlConfigDict | None

DEFAULT_QTL_TOOL_CONFIG: QtlToolConfigDict = {
    'n_bins': 100,
    'bed_fixed_fields': 2,
    'qtl_core_config': None
    }

@step(vtag='0.4 fix fix dofs mle')
def call_qtl(genotype_vcf: str, expression_bed: str, gt_covariates_file: str,
        gx_covariates_file: str, qtl_tool_config: QtlToolConfigDict | None = None):
    """
    Meta step to create the chunked qtl calling graph

    Args:
        genotype_vcf: path to the genotype in bgzipped vcf format
        expression_bed: path to expression
        gt_covariates_file: covariates to regress genotype on
        gx_covariates_file: covariates to regress expression on ; typically same
            as for genotype but may be changed if expression is pre-processed
            upstream of the qtl caller
        qtl_tool_config: dictionnary of options with keys:

            - **n_bins**: target number of bins into which the expressed genes
               will be split for parallelization. Final number of bins can
               differ, see :func:`gtex_qtl.qtl.split_genes`.
            - **bed_fixed_fields**: number of fields (leading columns) in the
                given `expression_bed` that are not samples. Must be at least
                two (chromosome and linear position of each feature must be the
                first two columns), but more fields can be present and will be
                added unmodified to the result files
            - **qtl_core_config**: dictionnary of further options to pass to
                :func:`gtex_qtl.qtl.call_qtls`

    """
    merged_qtl_tool_config = DEFAULT_QTL_TOOL_CONFIG  | (qtl_tool_config or {})

    expression_df = pd.read_table(expression_bed)

    windows = qtl.split_genes(expression_df, merged_qtl_tool_config['n_bins'])

    return merge_qtl([
        call_qtl_bin(genotype_vcf, expression_bed, gt_covariates_file,
            gx_covariates_file, window,
            merged_qtl_tool_config)
        for window in windows])

@step(vtag='0.3 fix fix dofs mle')
def call_qtl_bin(genotype_vcf: str, expression_bed: str,
        gt_covariates_file: str, gx_covariates_file: str,
        window: tuple[int, int], qtl_tool_config: QtlToolConfigDict):
    """
    Run the qtl caller on one chunk, and write down the results to disk
    """

    gt_covariates_df = pd.read_table(gt_covariates_file)
    if gx_covariates_file is None:
        gx_covariates_df = None
    else:
        gx_covariates_df = pd.read_table(gx_covariates_file)
    expression_df = pd.read_table(expression_bed)

    pairs, summary_perm, summary_ic = qtl.call_qtls(
            (expression_df, qtl_tool_config['bed_fixed_fields']),
            window, genotype_vcf,
            gt_covariates_df, gx_covariates_df,
            qtl_config=qtl_tool_config['qtl_core_config']
            )

    # pairs can be huge, so write it back and compress it
    pairs_path = galp.new_path()
    pairs.to_feather(pairs_path)

    # Gene level summaries are much lighter, so can use default handling
    return {
        'all_pairs_path': pairs_path,
        'egenes_perm': summary_perm,
        'egenes_ic': summary_ic
        }

@step(vtag='0.2 forward paths')
def merge_qtl(bins):
    """
    Concatenate results of all given bins
    """

#    FIXME: this is slow and temp-disabled for developping
#    pairs_path = galp.new_path()
#    with gzip.open(pairs_path, 'wt', encoding='utf8') as stream:
#        pd.read_feather(bins[0]['all_pairs_path']).to_csv(
#                stream, sep='\t', header=True
#                )
#        for _bin in bins[1:]:
#            pd.read_feather(_bin['all_pairs_path']).to_csv(
#                    stream, sep='\t', header=False
#                    )
#    pairs_path = None

    egenes_paths = {}
    for key in ('egenes_perm', 'egenes_ic'):
        path = galp.new_path()
        pd.concat([b[key] for b in bins], axis=0).to_csv(
                path, sep='\t', compression='gzip'
                )
        egenes_paths[key] = path

    return {
        'all_pairs': [b['all_pairs_path'] for b in bins],
        **egenes_paths
        }
