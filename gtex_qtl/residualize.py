"""
Pipeline fork demonstrating pre-residualization of expression
"""

from typing import TypedDict

import numpy as np
import numpy.typing as npt
import pandas as pd

import galp
from galp import step
import gemz_galp.models
import gemz.models

from . import stats

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

@step
def residualize(prepared_expression, covariates, model_spec):
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

@step
def _residualize_writeback(res_sg, samples, expr_meta):
    # Write back
    dst_path = galp.new_path()
    write_expression(dst_path, res_sg.T, samples, expr_meta)
    return dst_path

@step(items=2)
def split_covariates(combined_covariates_file):
    """
    Split the original combined covariate files into two sets.

    Returns:
        A tuple of paths to resp. the inferred and external covariates
    """
    cov_df = pd.read_table(combined_covariates_file)

    is_inf = cov_df['ID'].str.startswith('InferredCov')

    inf_df = cov_df[is_inf]
    ext_df = cov_df[~is_inf]

    inf_path = galp.new_path()
    inf_df.to_csv(inf_path, sep='\t', index=False)

    ext_path = galp.new_path()
    ext_df.to_csv(ext_path, sep='\t', index=False)

    return inf_path, ext_path

class CMKParamsDict(TypedDict):
    """
    Prameters optimized from CMK
    """
    groups_x: npt.NDArray
    compact_covariance_kk: npt.NDArray
    data_vars_x: npt.NDArray

@step(items=2)
def get_cmk_parameters(prepared_expression_path: str, covariates_path: str
        ) -> tuple[str, CMKParamsDict]:
    """
    Regress expression on covariates, optimize CMK parameters on residuals

    Returns the residuals (as a path) and the CMK parameters (as an inline dict
    of arrays).
    """
    # Read expression
    samples, expr_gs, expr_meta = read_expression(prepared_expression_path)

    # Read covariates
    cov_cs = read_covariates(samples, covariates_path)

    # Regress expression on covariates
    # Note: this will add an intercept before regressing
    res_gs = stats.regress(expr_gs, cov_cs)

    # Save residuals
    residuals_path = _residualize_writeback(res_gs.T, samples, expr_meta)

    # Run CMK/100
    cmk_fit = gemz.models.cmk.fit(res_gs.T, 100)

    if cmk_fit['aborted']:
        raise RuntimeError(f'CMK fit failed ({cmk_fit["aborted"]}): {cmk_fit["errors"]}')

    params: CMKParamsDict = {
            'groups_x': np.array(cmk_fit['data']['groups']),
            'compact_covariance_kk': np.array(cmk_fit['state']['compact_covariance']),
            'data_vars_x': np.array(cmk_fit['state']['data_vars'])
            }

    return residuals_path, params

def make_weight_function(params: CMKParamsDict):
    """
    Return a weight calculation function for each gene to use with the
    autoregression feature of the QTL caller
    """
    def get_weights(gene_index: int) -> tuple[float, npt.NDArray]:
        target_group = params['groups_x'][gene_index]
        # D1 is the predicted, D2 the predictor
        weights_x = np.take_along_axis(
                params['compact_covariance_kk'][target_group, :],
                params['groups_x'],
                axis=0)
        weights_x[gene_index] = 0
        # Data vars is a multiplier for both identity and groups, so we can just
        # drop it
        return 1.0, weights_x
    return get_weights
