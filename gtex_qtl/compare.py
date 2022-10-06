"""
Utils to compare fastqtl results
"""

import io
import base64

import pandas as pd
import datashader as ds
import datashader.transfer_functions as tf
import plotly.graph_objects as go

import galp

pbl = galp.StepSet()

ALL_PAIRS = 'fastqtl_workflow.fastqtl_nominal.allpairs'

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
def datashader_scatter(x, y):
    """
    Quick datashader scatter plot of two arrays.

    Arrays may be large but still need to fit in memory.
    """
    source = pd.DataFrame({
        'x': x,
        'y': y,
        })

    res = 4000 # 16 MP
    canvas = ds.Canvas(res, res)

    points = canvas.points(source, 'x', 'y')

    return tf.shade(points)

def plot_ds_scatter(ds_img):
    """
    Render a datashader scatter image in plotly
    """
    # Serialize image
    img_dict = {}
    with io.BytesIO() as buf:
        ds_img.to_pil().save(buf, format='png')
        img_dict['source'] = 'data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode('ascii')

    # Set position
    img_dict['x'] = ds_img['x'].data.min()
    img_dict['y'] = ds_img['y'].data.min()
    img_dict['sizex'] = ds_img['x'].data.max() - img_dict['x']
    img_dict['sizey'] = ds_img['y'].data.max() - img_dict['y']
    img_dict['yanchor'] = "bottom"
    img_dict['xref'] = 'x'
    img_dict['yref'] = 'y'
    img_dict['layer'] = 'below'


    range_x = img_dict['x'], img_dict['x'] + img_dict['sizex']
    range_y = img_dict['y'], img_dict['y'] + img_dict['sizey']

    # Generate figure
    fig = go.Figure(
        data=[
            go.Scatter(x=range_x, y=range_x, mode='lines',
                        line={'color': 'red'},
                        name='y = x', showlegend=True)
        ],
        layout={
            'plot_bgcolor': 'white',
            'height': 1000,
            'xaxis': { 'range': range_x,
                      'ticks': 'outside', 'linecolor': 'black',
                      'showgrid': False, 'zerolinecolor': 'black'},
            'yaxis': { 'range': range_y,
                      'scaleanchor': 'x', 'scaleratio': img_dict['sizex'] / img_dict['sizey'],
                     'ticks': 'outside', 'linecolor': 'black',
                      'showgrid': False, 'zerolinecolor': 'black'},
            'images': [img_dict]
        }
    )

    return fig

