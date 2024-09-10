"""
Collection of pipeline plots
"""

import os
import argparse
import numpy as np
import pandas as pd
import plotly.graph_objects as go

import galp
import gtex_qtl
from gtex_qtl import template

_EXCLUDES = {'linear'}

#pbl.bind(cmp_variable=
#        'qval'
        #'num_var'
        #'maf'
        #'pval_nominal'
        #'beta_shape1'
        #'beta_shape2'
        #'true_df'
        #'pval_beta'
#        )

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

def egenes_pval_cdf(all_egenes: dict[str, pd.DataFrame]) -> go.Figure:
    """
    CDF of adjusted p-values
    """
    return go.Figure([
        go.Scatter(
            x=np.sort(egenes['pval_beta']),
            y=np.arange(len(egenes)) + 1,
            mode='lines',
            name=pipeline,
            )
        for pipeline, egenes in sorted(all_egenes.items(),
            key=lambda p: p[1]['pval_beta'].median())
        if pipeline != 'linear'
        ], {
            #'title': 'CDF of gene-level adjusted p-values of best association, '
            # 'across all expressed genes',
            'xaxis': {
                'title': 'Multiple-comparison adjusted p-value of best '
                    'association',
                'type': 'log',
                'range': [-8.1, 0.1],
                },
            'yaxis': {
                'title': 'Number of genes (cumulative)',
                'exponentformat': 'none',
                'rangemode': 'tozero',
                },
            }
        )

def pairs_pval_cdf(all_pairs_pvals):
    """
    Quantiles of adjusted p-values of all associations
    """
    return go.Figure([
        go.Scatter(
            x=edges[1:],
            y=np.cumsum(counts),
            mode='lines',
            name=pipeline
            )
        for pipeline, (counts, edges) in sorted(
            all_pairs_pvals.items(),
            key=lambda p: -sum(p[1][0][:len(p[1][0]) // 2])
            )
        if pipeline not in _EXCLUDES
        ], {
            'title': 'Quantiles of gene-level adjusted p-values of all '
                'associations',
            'xaxis': {
                'title': 'P-value',
                'type': 'log',
                'range': [ -8.1, -0.9 ]
                },
            'yaxis': {
                'title': 'Number of associations (cumulative)',
                'exponentformat': 'SI',
                'rangemode': 'tozero',
                },
            }
        )

def vs_egenes_published_repro(all_egenes_perm: dict[str, pd.DataFrame]) -> go.Figure:
    """
    Plot number of eGenes
    """
    return vs_egenes(all_egenes_perm, (
        ('published', 'Published'),
        ('reproduction', 'Reproduction')
        ))
def vs_egenes_peer_cmk(all_egenes_ic: dict[str, pd.DataFrame]) -> go.Figure:
    """
    Plot number of eGenes
    """
    return vs_egenes(all_egenes_ic, (
        ('reimplementation', 'PEER'),
        ('cmk', 'cmk')
        ))

def vs_egenes(all_egenes: dict[str, pd.DataFrame], alts) -> go.Figure:
    """
    Plot number of eGenes
    """
    minexp = -5
    egenes = pd.DataFrame({
        name: [np.sum(all_egenes[name]['pval_beta'] < 10**k) for k in
               range(minexp, 0)]
        for name, _ in alts
    })

    # This needs to be tuned somewhat manually for best effect
    xmin, xmax = egenes[alts[0][0]].iloc[0] - 800, egenes[alts[1][0]].iloc[-1] + 400

    return go.Figure([
        go.Scatter({
            'x': egenes[alts[0][0]],
            'y': egenes[alts[1][0]],
            'text': [f'p < 10<sup>{k}</sup>' for k in range(minexp, 0)],
            'mode': 'lines+markers+text',
            'textposition': 'top left',
            'showlegend': False,
        }),
        go.Scatter({
            'x': [xmin, xmax],
            'y': [xmin, xmax],
            'mode': 'lines',
            'line': {'color': 'black', 'width': template.PX_RULE},
            'showlegend': False,
        })
    ], {
        'xaxis': {'range': [xmin, xmax],
            'title': f'Genes with significant eQTL ({alts[0][1]})'
            },
        'yaxis': {'range': [xmin, xmax],
        'title': f'Genes with significant eQTL ({alts[1][1]})'
        }
    })

def export_all_plots(store: str, output: str) -> None:
    """
    Run all the graph generating code and create svg files
    """
    os.makedirs(output, exist_ok=True)

    all_egenes_ic = galp.run(gtex_qtl.all_egenes_ic, store=store)
    all_egenes_perm = galp.run(gtex_qtl.all_egenes_perm, store=store)

    vs_egenes_peer_cmk(all_egenes_ic).write_image(
            os.path.join(output, 'egenes_peer_cmk.svg')
            )

    vs_egenes_published_repro(all_egenes_perm).write_image(
            os.path.join(output, 'egenes_pub_repro.svg')
            )

def main():
    """Entry point"""
    parser = argparse.ArgumentParser()
    parser.add_argument('store', help='Task result store')
    parser.add_argument('output', help='Task result store')
    args = parser.parse_args()

    export_all_plots(args.store, args.output)
