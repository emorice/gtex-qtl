"""
Wrapper around the qtl caller for pipelining
"""

import gzip
import galp
import pandas as pd

from . import qtl

pbl = galp.Block()


DEFAULT_QTL_TOOL_CONFIG = {
    'n_bins': 100,
    'bed_fixed_fields': 2,
    'qtl_core_config': None
    }

@pbl.step(vtag='0.2 merge_qtl-0.2')
def call_qtl(genotype_vcf, expression_bed, gt_covariates_file,
        gx_covariates_file, qtl_tool_config=None):
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
    qtl_tool_config = dict(DEFAULT_QTL_TOOL_CONFIG, **qtl_tool_config)

    expression_df = pd.read_table(expression_bed)

    windows = qtl.split_genes(expression_df, qtl_tool_config['n_bins'])

    return merge_qtl([
        call_qtl_bin(genotype_vcf, expression_bed, gt_covariates_file,
            gx_covariates_file, window,
            qtl_tool_config)
        for window in windows])

@pbl.step
def call_qtl_bin(genotype_vcf, expression_bed, gt_covariates_file,
        gx_covariates_file, window, qtl_tool_config, _galp):
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
    pairs_path = _galp.new_path()
    pairs.to_feather(pairs_path)

    # Gene level summaries are much lighter, so can use default handling
    return {
        'all_pairs_path': pairs_path,
        'egenes_perm': summary_perm,
        'egenes_ic': summary_ic
        }

@pbl.step(vtag='0.2 forward paths')
def merge_qtl(bins, _galp):
    """
    Concatenate results of all given bins
    """

#    FIXME: this is slow and temp-disabled for developping
#    pairs_path = _galp.new_path()
#    with gzip.open(pairs_path, 'wt', encoding='utf8') as stream:
#        pd.read_feather(bins[0]['all_pairs_path']).to_csv(
#                stream, sep='\t', header=True
#                )
#        for _bin in bins[1:]:
#            pd.read_feather(_bin['all_pairs_path']).to_csv(
#                    stream, sep='\t', header=False
#                    )
    pairs_path = None

    egenes_paths = {}
    for key in ('egenes_perm', 'egenes_ic'):
        path = _galp.new_path()
        pd.concat([b[key] for b in bins], axis=0).to_csv(
                path, sep='\t', compression='gzip'
                )
        egenes_paths[key] = path

    return {
        'all_pairs': [b['all_pairs_path'] for b in bins],
        **egenes_paths
        }
