"""
Utils to compare fastqtl results
"""

import io
import base64

import numpy as np
import pandas as pd
import datashader as ds
import datashader.transfer_functions as tf
import plotly.graph_objects as go

import galp

pbl = galp.Block()

ALL_PAIRS = 'fastqtl_workflow.fastqtl_nominal.allpairs'
EGENES = 'fastqtl_workflow.fastqtl_postprocess.genes_annotated'
SIGNIF_PAIRS = 'fastqtl_workflow.fastqtl_postprocess.signifpairs'

@pbl.step
def all_pvals(results):
    """
    Extract all nominal p-values as an array
    """
    return (
        pd.read_table(
                results[ALL_PAIRS],
                usecols=['pval_nominal'],
                dtype='float32'
            )
        ['pval_nominal']
        .to_numpy()
    )

@pbl.step
def histogram(values, log=False, nbins=1000):
    """
    Pre-compute an histogram.

    Silently discards infs/nans.
    """
    trans = np.log10 if log else (lambda x: x)

    values = trans(values)

    values = values[np.isfinite(values)]

    hist = np.histogram(
            values,
            bins=nbins,
            )

    return {
        'counts': hist[0],
        'edges': hist[1],
        'log': log
        }

@pbl.step
def quantiles(values, num=1000):
    """
    Precompute quantiles
    """
    probas = np.linspace(0., 1, num+1)
    quantiles = np.quantile(values, probas)
    return probas, quantiles

def plot_histogram(hist):
    """
    Render a pre-computed histogram
    """

    counts = hist['counts']
    if hist['log']:
        edges = np.exp(hist['edges'])
    else:
        edges = hist['edges']
    middles = .5 * (edges[1:] + edges[:-1])

    return go.Figure(
            data=[
                go.Bar(
                    y=counts,
                    x=middles,
                    width=(edges[1:] - edges[:-1]) * 0.9,
                    marker={'line': {'width': 0}}
                    )
                ],
            layout={
                'xaxis': {
                    'type': 'log' if hist['log'] else 'linear',
                    },
                'width': 1200,
                }
            )

@pbl.step(vtag='0.4: trend')
def datashader_scatter(x, y, log=False):
    """
    Quick datashader scatter plot of two arrays.

    Arrays may be large but still need to fit in memory.
    """
    trans = np.log10 if log else (lambda x: x)

    source = (pd.DataFrame({
        'x': trans(x),
        'y': trans(y),
        })
        .replace(-np.inf, np.nan)
        .dropna()
        )

    n_bins = 1000
    trend = (
        source
        .sort_values('x')
        .assign(cdf=
            np.floor(n_bins * np.arange(len(source)) / len(source))
            / n_bins
            )
        .groupby('cdf')
        .mean()
        )

    res = 4000 # 16 MP
    canvas = ds.Canvas(res, res)

    points = canvas.points(source, 'x', 'y')

    return tf.shade(points), trend

def plot_ds_scatter(ds_image_trend):
    """
    Render a datashader scatter image in plotly
    """
    ds_img, trend = ds_image_trend
    # Serialize image
    img_dict = {}
    with io.BytesIO() as buf:
        ds_img.to_pil().save(buf, format='png')
        img_dict['source'] = (
                'data:image/png;base64,' +
                base64.b64encode(buf.getvalue()).decode('ascii')
                )

    # Set position
    img_dict['x'] = ds_img['x'].data.min()
    img_dict['y'] = ds_img['y'].data.min()
    img_dict['sizex'] = ds_img['x'].data.max() - img_dict['x']
    img_dict['sizey'] = ds_img['y'].data.max() - img_dict['y']
    img_dict['yanchor'] = "bottom"
    img_dict['xref'] = 'x'
    img_dict['yref'] = 'y'
    img_dict['layer'] = 'below'


    margin_ratio = 0.05
    ranges = { k: (
            img_dict[k] - margin_ratio * img_dict[f'size{k}'],
            img_dict[k] + (1. + margin_ratio) * img_dict[f'size{k}']
            ) for k in 'xy'
            }

    # Generate figure
    fig = go.Figure(
        data=[
            go.Scatter(x=ranges['x'], y=ranges['x'], mode='lines',
                line={'color': 'red', 'width': 1},
                        name='y = x', showlegend=True),
            go.Scatter(x=trend['x'], y=trend['y'], mode='lines',
                        hovertext=[
                            f'{100 * l:.6g} - {100 *u:.6g} %'
                            for l, u in
                            zip(trend.index, [*trend.index[1:], 1.0])
                            ],
                        line={'color': 'green'},
                        name='mean', showlegend=True)
        ],
        layout={
            'plot_bgcolor': 'white',
            'height': 1000,
            'xaxis': { 'range': ranges['x'],
                      'ticks': 'outside', 'linecolor': 'black',
                      'showgrid': False, 'zerolinecolor': 'black'},
            'yaxis': { 'range': ranges['y'],
                      'scaleanchor': 'x', 'scaleratio': img_dict['sizex'] / img_dict['sizey'],
                     'ticks': 'outside', 'linecolor': 'black',
                      'showgrid': False, 'zerolinecolor': 'black'},
            'images': [img_dict]
        }
    )

    return fig

@pbl.step(vtag='0.2 pairs')
def count_egenes(results):
    """
    Count genes with global qval <= 0.05, and significant eqtls.
    """
    return (
        len(
            pd.read_table(
                results[EGENES]
                )
            .query("qval <= 0.05")
            ),
        len(
            pd.read_table(
                results[SIGNIF_PAIRS]
                )
            )
        )

@pbl.step
def all_egenes(egenes, all_qtls):
    """
    Extract gene summaries computed by permutation
    """

    egene_files = dict(egenes, **{
        pipeline: (
            results[EGENES] if EGENES in results
            else results['egenes_perm']
            )
        for pipeline, results in all_qtls.items()
        })

    return {
        pipeline: pd.read_table(file, compression='gzip')
        for pipeline, file in egene_files.items()
        }
