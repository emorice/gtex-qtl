"""
Wrapper around the qtl caller for pipelining
"""

from typing import TypedDict

import pandas as pd
import galp
from galp import step

from . import qtl, residualize

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

@step(vtag='0.6 accept missing gx covs')
def call_qtl(genotype_vcf: str, expression_bed: str, gt_covariates_file: str,
        gx_covariates_file: str | None, qtl_tool_config: QtlToolConfigDict | None =
        None, autoregress: bool | dict = False):
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
        autoregress: whether to include leave-one-gene-out regression in the
            controls, for both genotype and expression. Can be a dict with a key
            indicating the precise method and value the auxiliary data for said
            method.

    """
    merged_qtl_tool_config = DEFAULT_QTL_TOOL_CONFIG  | (qtl_tool_config or {})

    expression_df = pd.read_table(expression_bed)

    windows = qtl.split_genes(expression_df, merged_qtl_tool_config['n_bins'])

    return merge_qtl([
        call_qtl_bin(genotype_vcf, expression_bed, gt_covariates_file,
            gx_covariates_file, window,
            merged_qtl_tool_config, autoregress)
        for window in windows])

@step(vtag='0.5 accept missing gx covs')
def call_qtl_bin(genotype_vcf: str, expression_bed: str,
        gt_covariates_file: str, gx_covariates_file: str | None,
        window: tuple[int, int], qtl_tool_config: QtlToolConfigDict,
        autoregress: bool):
    """
    Run the qtl caller on one chunk, and write down the results to disk
    """

    gt_covariates_df = pd.read_table(gt_covariates_file)
    if gx_covariates_file is None:
        gx_covariates_df = None
    else:
        gx_covariates_df = pd.read_table(gx_covariates_file)
    expression_df = pd.read_table(expression_bed)

    gt_regressors: list[qtl.Regressor] = [
            {'method': 'external', 'data': gt_covariates_df}]

    gx_regressors: list[qtl.Regressor] = []
    if gx_covariates_df is not None:
        gx_regressors.append(
            {'method': 'external', 'data': gx_covariates_df}
            )
    if autoregress:
        if autoregress is True:
            auto_data = None
        elif isinstance(autoregress, dict) and 'cmk_params' in autoregress:
            auto_data = residualize.make_weight_function(
                autoregress['cmk_params']
                )
        else:
            raise ValueError(f'Bad autoregress strategy: {autoregress!r}')

        gt_regressors.append({'method': 'auto', 'data': auto_data})
        gx_regressors.append({'method': 'auto', 'data': auto_data})
    config = (qtl_tool_config['qtl_core_config'] or {}).copy()
    config['genotype_regressors'] = gt_regressors
    config['expression_regressors'] = gx_regressors

    pairs, summary_perm, summary_ic = qtl.call_qtls(
            (expression_df, qtl_tool_config['bed_fixed_fields']),
            window, genotype_vcf, qtl_config=config
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
