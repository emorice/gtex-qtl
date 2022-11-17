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

# 4) Run fastqtl
# ==============

run_fastqtl = functools.partial(fastqtl.run_fastqtl,
        indexed_vcf=preprocess.filtered_indexed_vcf,
        gene_model=public_input_files['gene_model'],
        )

reference_fastqtl = run_fastqtl(
        prepared_expression,
        covariates_file=combined_covariates_file
        )

# 5) Get published results and compare
# ====================================

plots.pbl.bind(published_tissue_egenes=downloads.published_tissue_egenes())

plots.pbl.bind(computed_tissue_egenes=
        reference_fastqtl['fastqtl_workflow.fastqtl_permutations_merge.genes']
        )

# 6) Pre-residualization alternatives
# ===================================

# 6.1) Simplest pre-residualization
# ---------------------------------

pre_covariates_file, post_covariate_file = (
    residualize.split_covariates(combined_covariates_file)
    )

residualized_expression = preprocess.tabix(
        residualize.residualize(
            prepared_expression[0],
            pre_covariates_file)
        )

residualized_fastqtl = run_fastqtl(residualized_expression,
            covariates_file=post_covariate_file)

# 6.2) Other models
# ------------------

model_specs = {
        'linear': {'model': 'linear'},
        'peer': {'model': 'peer', 'n_factors': 60},
        'cmk': {'model': 'cmk', 'n_groups': 100},
        }

blind_expressions = {
        shorthand: preprocess.tabix(
            residualize.residualize_blind(
                prepared_expression[0],
                model_spec
                )
            )
        for shorthand, model_spec in model_specs.items()
        }

blind_fastqtl = {
        shorthand: run_fastqtl(
            blind_expression,
            covariates_file=post_covariate_file)
        for shorthand, blind_expression in blind_expressions.items()
        }

# 6.3) Gather runs
# ----------------

all_fastqtl = {
        'reproduction': fastqtl,
        'pre-residualization': residualized_fastqtl,
        **blind_fastqtl
        }

# 6.3) Plots
# ----------

pbl = galp.Block()

pbl.bind(res_vs_orig_raster=compare.datashader_scatter(
        compare.all_pvals(fastqtl),
        compare.all_pvals(residualized_fastqtl),
        log=False
        ))

@pbl.view
def residualized_pvals_plot(res_vs_orig_raster):
    """
    Residualized vs original p-values for all pairs
    """
    fig = compare.plot_ds_scatter(res_vs_orig_raster)
    fig.update_layout({
        'title':
        'Comparison of nominal p-values for all tested gene-variant pairs',
        'xaxis.title': 'Reproduced p-value',
        'yaxis.title': 'Pre-residualized p-value',
        })
    return fig


pbl.bind(blind_vs_res_raster=compare.datashader_scatter(
        compare.all_pvals(residualized_fastqtl),
        compare.all_pvals(blind_fastqtl['linear']),
        log=True
        ))

@pbl.view
def blind_pvals_plot(blind_vs_res_raster):
    """
    Blind linear vs residualized p-values for all pairs
    """
    fig = compare.plot_ds_scatter(blind_vs_res_raster)
    fig.update_layout({
        'title':
        'Comparison of nominal p-values for all tested gene-variant pairs',
        'xaxis.title': 'Pre-residualized p-value',
        'yaxis.title': 'Blind linear p-value',
        })
    return fig

pbl.bind(pvals_hist=compare.histogram(
    compare.all_pvals(residualized_fastqtl),
    ))

@pbl.view
def pvals_histogram_plot(pvals_hist):
    """
    Histogram of p-values for all pairs
    """
    fig = compare.plot_histogram(pvals_hist)
    fig.update_layout({
        'title':
        'Distribution of raw p-values for all tested gene-variant pairs',
        'xaxis.title': 'p-value',
        'yaxis.title': 'Number of gene-variant pairs'
        })
    return fig

pbl.bind(pvals_quantiles_ref=compare.quantiles(
    compare.all_pvals(all_fastqtl['pre-residualization'])
    ))
pbl.bind(pvals_quantiles_alt=compare.quantiles(
    compare.all_pvals(all_fastqtl['cmk'])
    ))

pbl.bind(qq_log=True)
pbl.bind(qq_relative=True)

@pbl.view
def pvals_qq_plot(pvals_quantiles_alt, pvals_quantiles_ref=None, qq_log=False,
        qq_relative=False):
    """
    Uniform qq plot
    """

    probas_alt, quantiles_alt = pvals_quantiles_alt

    if pvals_quantiles_ref is None:
        probas_ref, quantiles_ref = probas_alt, probas_alt
    else:
        probas_ref, quantiles_ref = pvals_quantiles_ref

    if not np.allclose(probas_ref, probas_alt):
        raise ValueError('Need matching probas')

    neutral = np.array(quantiles_ref)
    if qq_relative:
        quantiles_alt /= quantiles_ref
        neutral /= quantiles_ref

    return go.Figure(
        data=[
            go.Scatter(
                x=quantiles_ref,
                y=quantiles_alt,
                mode='lines',
                name='quantiles',
                hovertext=[f'{100 * p:.6g} %' for p in probas_ref],
                yaxis='y2'
            ),
            go.Scatter(
                x=quantiles_ref,
                y=neutral,
                mode='lines',
                name='y = x',
                yaxis='y2'
                ),
            go.Scatter(
                x=quantiles_ref,
                y=probas_ref,
                mode='lines',
                name='Reference CDF',
                yaxis='y1'
                )
            ],
        layout={
            'width': 1000,
            'height': 1250,
            'xaxis': {
                'title': 'Reference p-values quantiles',
                'type': 'log' if qq_log else 'linear',
                'exponentformat': 'power',
                'linecolor': 'black',
                'ticks': 'outside',
                'position': 0.22,
                },
            'yaxis2': {
                'title': 'Alternative p-values quantiles',
                'type': 'log' if qq_log and not qq_relative else 'linear',
                'exponentformat': 'power',
                'linecolor': 'black',
                'ticks': 'outside',
                'domain': [0.22, 1.0]
                },
            'yaxis1': {
                'title': 'CDF',
                'domain': [0.0, 0.16],
                'type': 'log',
                'ticks': 'outside',
                'linecolor': 'black',
                },
            }
        )

pbl.bind(egene_counts={
    exp_name: compare.count_egenes(results)
    for exp_name, results in all_fastqtl.items()
    })

@pbl.view
def egenes_plot(egene_counts):
    """
    Bar chart of number of "eGenes"
    """

    names, counts = zip(*egene_counts.items())
    egenes, pairs = zip(*counts)

    return go.Figure(
            data=[
                go.Bar(
                    x=egenes,
                    y=names,
                    orientation='h',
                    name='eGenes',
                    ),
                go.Bar(
                    x=pairs,
                    y=names,
                    orientation='h',
                    xaxis='x2',
                    name='gene-variant pairs',
                    )
            ],
            layout={
                'title': 'Number of "eGenes" (genes with one significant eQTL '
                    'at FDR ≤ 0.05) and gene-variant pairs per method',
                    'yaxis': {'title': 'Pipeline run code name'},
                    'xaxis': {'title': 'Number of genes',
                        'rangemode': 'tozero',
                        'exponentformat': 'none',
                        'domain': [0., 0.49]},
                    'xaxis2': {'title': 'Number of significant variant-genes pairs',
                        'rangemode': 'tozero',
                        'exponentformat': 'none',
                        'domain': [0.51, 1.0]},
                }
            )

# END
# ===

#plots = [ residualized_pvals_plot, blind_pvals_plot ]

default_target = plots.qvalues_cmp # reference_fastqtl # egenes_plot # all_fastqtl
