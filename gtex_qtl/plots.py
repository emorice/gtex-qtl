"""
Collection of pipeline plots
"""

import numpy as np
import galp
import pandas as pd
import plotly.graph_objects as go

pbl = galp.Block()


pbl.bind(cmp_variable=
        'qval'
        #'num_var'
        #'maf'
        #'pval_nominal'
        #'beta_shape1'
        #'beta_shape2'
        #'true_df'
        #'pval_beta'
        )

@pbl.view
def qvalues_cmp(published_tissue_egenes, computed_tissue_egenes,
    cmp_variable='qval'):
    """
    Scatter plot of egenes variables
    """
    merged = (
        pd.read_table(published_tissue_egenes)
        .merge(
            pd.read_table(computed_tissue_egenes),
            on='gene_id',
            suffixes=['_published', '_computed']
            )
        )

    pub = merged[cmp_variable + '_published']
    alt = merged[cmp_variable + '_computed']

    return go.Figure(
        data=[
            go.Scattergl(
                x=pub,
                y=pub,
                mode='lines',
                name='y = x',
                line={'width': 2, 'color': 'grey'},
                yaxis='y2'
                ),
            go.Scattergl(
                x=pub,
                y=alt,
                mode='markers',
                marker={'size': 3, 'color': 'blue'},
                name='recomputed',
                yaxis='y2'
                ),
            go.Histogram(
                x=pub,
                yaxis='y1'
                )
            ],
        layout={
            'title': 'Comparison of published and reproduced gene statistics',
            'xaxis': {'type': 'linear', 'title': f'Gene-level {cmp_variable} (published)'},
            'yaxis2': {
                'type': 'linear',
                'title': f'Gene-level {cmp_variable} (recomputed)',
                'domain': [0.21, 1.0]
                },
            'yaxis1': {
                'domain': [0.0, 0.19]
                },
            'height': 1000,
            'width': 1000,
            }
        )

@pbl.view
def egenes_pval_cdf(all_egenes):
    """
    CDF of adjusted p-values
    """
    return go.Figure([
        go.Scatter(
            x=np.sort(egenes['pval_beta']),
            y=np.arange(len(egenes)) + 1,
            mode='lines',
            name=pipeline
            )
        for pipeline, egenes in all_egenes.items()
        ], {
            'title': 'CDF of gene-level adjusted pvalues of best association, '
                'across all expressed genes',
            'xaxis': {
                'title': 'Multiple-comparison-adjusted pvalue of best '
                    'association',
                'type': 'log',
                'exponentformat': 'power',
                },
            'yaxis': {
                'title': 'Number of genes (cumulative)'
                },
            'height': 1000,
            }
        )
