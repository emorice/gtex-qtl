"""
GTEx QTL calling pipeline

This file defines the reference reproduction pipeline. Sister modules contain
pipeline variants.
"""

import functools

from . import (
        downloads,
        preprocess,
        peer,
        fastqtl,
        residualize, compare,
        qtl_tool,
        )

# pylint: disable=redefined-outer-name


def main() -> tuple[dict, dict]:
    """
    Main pipeline.

    Returns two tasks yielding dictionaries of dataframes, corresponding to the gene-level
    summaries with permutation or inter-chromosome control, respectively, for
    all pipeline variants.
    """

    # 0.1) Extract additional (non-PEER) covariates and genotyped subject list
    # ========================================================================

    additional_covariates, genotyped_subject_ids = preprocess.get_additional_covariates()

    # 0.2) Download public expression files and gene model
    # ====================================================

    public_input_files = downloads.input_files

    for input_file in 'wb_tpm', 'wb_counts':
        public_input_files[input_file] = preprocess.gct_filter_columns(
            public_input_files[input_file],
            genotyped_subject_ids
            )

    # 1) Prepare expression (filter and normalize)
    # ============================================

    filtered_indexed_vcf = preprocess.filter_index_vcf()

    prepared_expression = preprocess.prepare_expression(
            public_input_files['wb_tpm'],
            public_input_files['wb_counts'],
            public_input_files['gene_model'],
            filtered_indexed_vcf
            )

    # 2) PEER factors and 3) combine covariates
    # =========================================

    combined_covariates_file = peer.peer_and_combine(
            prepared_expression[0],
            preprocess.gpcs_covariates(),
            additional_covariates)

    # 4) Run fastqtl and the integrated qtl caller
    # ============================================

    run_fastqtl = functools.partial(fastqtl.run_fastqtl,
            indexed_vcf=filtered_indexed_vcf,
            gene_model=public_input_files['gene_model'],
            )

    reference_fastqtl = run_fastqtl(
            prepared_expression,
            covariates_file=combined_covariates_file
            )

    # Fix common options for all pipelines run that follow
    run_qtl = functools.partial(qtl_tool.call_qtl,
            genotype_vcf=preprocess.index_vcf()[0],
            qtl_tool_config={
                'bed_fixed_fields': 4,
                'qtl_core_config': {
                    # Use small number of permutations for testing only
                    'num_null_genes': 1000
                    }
                }
            )

    reference_qtl = run_qtl(
            expression_bed=prepared_expression[0],
            gt_covariates_file=combined_covariates_file,
            gx_covariates_file=combined_covariates_file,
            )

    # 6) Alternative control strategies
    # ===================================

    # Inferred: PEERs, External: Genotypic PCs, Sequencing, Sex
    _inf_covariates_file, ext_covariates_file = (
        residualize.split_covariates(combined_covariates_file)
        )

    # 6.1) Other models
    # ------------------

    model_specs = {
            #'linear': {'model': 'linear'},
            'peer': {'model': 'peer', 'n_factors': 60},
            'cmk': {'model': 'cmk', 'n_groups': 100},
            }

    # Regress first on external covariates, then with gene covariance model
    residualized_expressions = {
            shorthand: residualize.residualize(
                prepared_expression[0],
                ext_covariates_file,
                model_spec
                )
            for shorthand, model_spec in model_specs.items()
            }

    # Regress only on external covariates, then fit model to obtain parameters
    ext_residualized_expr, cmk_params = residualize.get_cmk_parameters(
            prepared_expression[0],
            ext_covariates_file,
            )

    # Pipeline variants with pre-residualisation
    alt_qtl = {
            shorthand: run_qtl(
                expression_bed=expression,
                gt_covariates_file=combined_covariates_file, # Regress GT on Ext+PEER
                gx_covariates_file=None, # OTOH assume GX is already regressed
                )
            for shorthand, expression in residualized_expressions.items()
            }

    # Pipeline variant with autoregression inside caller
    alt_qtl['auto'] = run_qtl(
                expression_bed=prepared_expression[0],
                gt_covariates_file=ext_covariates_file, # regress gt in Ext only
                gx_covariates_file=ext_covariates_file, # regress gx in Ext only
                autoregress=True, # regress both on all genes but target
                )

    # Pipeline variants with CMK autoregression inside caller
    alt_qtl['auto_cmk'] = run_qtl(
                expression_bed=ext_residualized_expr, # Already regressed expression
                gt_covariates_file=ext_covariates_file, # regress gt in Ext only
                gx_covariates_file=None, # Gx is already regressed
                autoregress={'cmk_params': cmk_params}, # regress both with CMK weights
                )

    # 6.3) Gather runs
    # ----------------

    all_qtl = {
            'reproduction': reference_fastqtl,
            'reimplementation': reference_qtl,
            **alt_qtl
            }

    all_egenes_perm = {
        'published': downloads.get_published_tissue_egenes('Whole_Blood'),
        'reproduction': all_qtl['reproduction'][compare.EGENES],
        'reimplementation (perm)': all_qtl['reimplementation']['egenes_perm'],
        'reimplementation (interchrom)': all_qtl['reimplementation']['egenes_ic'],
        'autoregression': all_qtl['auto']['egenes_perm'],
        'autoregression (cmk)': all_qtl['auto_cmk']['egenes_perm'],
        }

    # All runs with comparable inter-chrom statistics
    qtls_ic = {'reimplementation': reference_qtl, **alt_qtl}

    # Gene-level summaries for all inter-chrom runs
    all_egenes_ic = {
        k: qtl['egenes_ic'] for k, qtl in qtls_ic.items()
        }


    return all_egenes_perm, all_egenes_ic

default_target = main
"""
For galp command line interface
"""
