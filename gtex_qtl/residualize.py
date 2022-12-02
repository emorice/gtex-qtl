"""
Pipeline fork demonstrating pre-residualization of expression
"""

import numpy as np
import pandas as pd

import galp
import gemz_galp.models

from . import stats

pbl = galp.Block()

def read_expression(path):
    """
    Load expression file
    """
    expr_df_gs = pd.read_table(path)

    samples = [ col for col in expr_df_gs if col.startswith('GTEX-') ]
    expr_gs = expr_df_gs[samples].to_numpy()
    expr_meta = expr_df_gs.drop(samples, axis=1)

    return samples, expr_gs, expr_meta

def read_covariates(samples, *paths):
    """
    Load covariate files
    """
    cov_cs_list = []

    for path in paths:
        cov_df_cs = pd.read_table(path)

        # Format as matching arrays
        cov_cs_list.append(cov_df_cs[samples].to_numpy())

    return np.vstack(cov_cs_list)

def write_expression(path, expr_gs, samples, expr_meta):
    """
    Write back expression
    """
    # Reformat as dataframe
    expr_df_gs = pd.DataFrame(expr_gs, columns=samples)

    # Reinclude other columns
    expr_df_gs = pd.concat((
        expr_meta,
        expr_df_gs
        ), axis=1)

    # Write back as uncompressed tsv
    expr_df_gs.to_csv(path, sep='\t', index=False)

@pbl.step
def residualize(prepared_expression, covariates, model_spec, _galp):
    """
    Regress expression first on covariates, then on itself with the specified
    model
    """
    # Read expression
    samples, expr_gs, expr_meta = read_expression(prepared_expression)

    # Read covariates
    cov_cs = read_covariates(samples, covariates)

    # Regress expression on covariates
    # Note: this will add an intercept before regressing
    res_gs = stats.regress(expr_gs, cov_cs)

    # Build parallel cross-regression graph
    return _residualize_writeback(
        gemz_galp.models.cv_residualize(model_spec, res_gs.T),
        samples, expr_meta)

@pbl.step
def _residualize_writeback(res_sg, samples, expr_meta, _galp):
    # Write back
    dst_path = _galp.new_path()
    write_expression(dst_path, res_sg.T, samples, expr_meta)
    return dst_path

@pbl.step(items=2)
def split_covariates(combined_covariates_file, _galp):
    """
    Split the original combined covariate files into two sets.

    Returns:
        A tuple of paths to resp. the inferred and external covariates
    """
    cov_df = pd.read_table(combined_covariates_file)

    is_inf = cov_df['ID'].str.startswith('InferredCov')

    inf_df = cov_df[is_inf]
    ext_df = cov_df[~is_inf]

    inf_path = _galp.new_path()
    inf_df.to_csv(inf_path, sep='\t', index=False)

    ext_path = _galp.new_path()
    ext_df.to_csv(ext_path, sep='\t', index=False)

    return inf_path, ext_path
