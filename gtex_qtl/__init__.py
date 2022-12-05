"""
GTEx QTL calling pipeline

This file defines the reference reproduction pipeline. Sister modules contain
pipeline variants.
"""

import os
import re
import gzip
import contextlib
import urllib.request
import shutil
import subprocess
import functools

import numpy as np
import pandas as pd
import plotly.graph_objects as go

import galp
import wdl_galp

from . import (
        downloads,
        preprocess,
        peer,
        fastqtl,
        residualize, compare,
        plots,
        qtl_tool,
        )

# pylint: disable=redefined-outer-name


# 0.1) Extract additional (non-PEER) covariates and genotyped subject list
# ========================================================================

additional_covariates, genotyped_subject_ids = preprocess.additional_covariates()

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

prepared_expression = preprocess.prepare_expression(
        public_input_files['wb_tpm'],
        public_input_files['wb_counts'],
        public_input_files['gene_model']
        )

# 2) PEER factors and 3) combine covariates
# =========================================

combined_covariates_file = peer.peer_and_combine(
        prepared_expression[0],
        preprocess.gpcs_covariates,
        additional_covariates)

# 4) Run fastqtl and the integrated qtl caller
# ============================================

run_fastqtl = functools.partial(fastqtl.run_fastqtl,
        indexed_vcf=preprocess.filtered_indexed_vcf,
        gene_model=public_input_files['gene_model'],
        )

reference_fastqtl = run_fastqtl(
        prepared_expression,
        covariates_file=combined_covariates_file
        )

run_qtl = functools.partial(qtl_tool.call_qtl,
        genotype_vcf=preprocess.indexed_vcf[0],
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

# 5) Get published results and compare
# ====================================

plots.pbl.bind(published_tissue_egenes=downloads.published_tissue_egenes())

plots.pbl.bind(computed_tissue_egenes=
        reference_fastqtl['fastqtl_workflow.fastqtl_permutations_merge.genes']
        )

# Models and plots below this point are undergoing rewriting

# 6) Residualization alternatives
# ===================================

# Inferred: PEERs, External: GPCS, Sequencing, Sex
inf_covariates_file, ext_covariates_file = (
    residualize.split_covariates(combined_covariates_file)
    )

# 6.1) Other models
# ------------------

model_specs = {
        'linear': {'model': 'linear'},
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

alt_qtl = {
        shorthand: run_qtl(
            expression_bed=expression,
            gt_covariates_file=combined_covariates_file, # Regress GT on Ext+PEER
            gx_covariates_file=None, # OTOH assume GX is already regressed
            )
        for shorthand, expression in residualized_expressions.items()
        }

# alt_qtl['linear_nogtpeer100'] = qtl_tool.call_qtl(
#         genotype_vcf=preprocess.indexed_vcf[0],
#         qtl_tool_config={
#             'bed_fixed_fields': 4,
#             'qtl_core_config': {
#                 # Use small number of permutations for testing only
#                 'num_null_genes': 100
#                 }
#             },
#         expression_bed=residualized_expressions['linear'],
#         gt_covariates_file=ext_covariates_file, # Regress GT on Ext only
#         gx_covariates_file=None, # OTOH assume GX is already regressed
#         )

# 6.3) Gather runs
# ----------------

all_qtl = {
        'reproduction': reference_fastqtl,
        'reimplementation': reference_qtl,
        **alt_qtl
        }

all_egenes_perm = compare.all_egenes({
    'published': downloads.published_tissue_egenes(),
    'reproduction': all_qtl['reproduction'][compare.EGENES],
    'reimplementation (perm)': all_qtl['reimplementation']['egenes_perm'],
    'reimplementation (interchrom)': all_qtl['reimplementation']['egenes_ic'],
    })

qtls_ic = {'reimplementation': reference_qtl, **alt_qtl}
"""
All runs with comparable inter-chrom statistics
"""

all_egenes_ic = compare.all_egenes({
    k: qtl['egenes_ic'] for k, qtl in qtls_ic.items()
    })
"""
Gene-level summaries for all inter-chrom runs
"""

plots.pbl.bind(all_egenes=all_egenes_ic)

all_pval_hist_ic = {
        ppl: compare.all_pairs_adjusted_hist(qtl)
        for ppl, qtl in qtls_ic.items()
        }

plots.pbl.bind(all_pairs_pvals=all_pval_hist_ic)

# END
# ===

default_target = all_pval_hist_ic
