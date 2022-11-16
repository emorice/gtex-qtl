"""
Pipeline fork demonstrating pre-residualization of expression
"""

import numpy as np
import pandas as pd

import galp
import gemz.models

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

def regress(expr_gs, cov_cs):
    """
    Regress expression against covariates and return residuals
    """
    # Add the intercept pseudo-covariate
    cov_cs = np.vstack((np.ones(expr_gs.shape[-1]), cov_cs))

    # Compute
    res_gs = (expr_gs - (
        (expr_gs @ cov_cs.T) @ np.linalg.inv(cov_cs @ cov_cs.T) @ cov_cs
        ))

    return res_gs

@pbl.step(vtag='0.2: intercept')
def residualize(prepared_expression, combined_covariates, _galp):
    """
    Perform the residualization as a separate step and write back the result as
    a new expression file.

    Here residualization is taken as meaning regressing the prepared expression
    linearily on the covariates and taking the residuals of this regression.
    """

    # Shorthands: g - genes, s - samples, c - covariates

    # Read expression
    samples, expr_gs, expr_meta = read_expression(prepared_expression)

    # Read covariates
    cov_cs = read_covariates(samples, combined_covariates)

    res_gs = regress(expr_gs, cov_cs)

    # Write back as uncompressed tsv
    dst_path = _galp.new_path() + '.bed'
    write_expression(dst_path, res_gs, samples, expr_meta)
    return dst_path

@pbl.step(vtag='0.2: intercept')
def residualize_blind(prepared_expression, model_spec, _galp):
    """
    Regress on given covariates, then linearily on other genes
    """

    # Read expression
    samples, expr_gs, expr_meta = read_expression(prepared_expression)

    # Perform main blind regression
    #model_spec = {'model': 'linear'}
    ## small dimension first
    res_gs = gemz.models.cv_residualize(model_spec, expr_gs.T).T

    # Write back
    dst_path = _galp.new_path() + '.bed'
    write_expression(dst_path, res_gs, samples, expr_meta)
    return dst_path

@pbl.step(items=2)
def split_covariates(combined_covariates_file, _galp):
    """
    Split the original combined covariate files into two sets.

    Returns:
        A tuple of paths to resp. the pre and post covariate files.

    The first set are covariates to use in pre-resisdualization, the second the
    set to pass to the QTL caller.
    """
    cov_df = pd.read_table(combined_covariates_file)

    is_pre = cov_df['ID'].str.startswith('InferredCov')

    pre_df = cov_df[is_pre]
    post_df = cov_df[~is_pre]

    pre_path = _galp.new_path()
    pre_df.to_csv(pre_path, sep='\t', index=False)

    post_path = _galp.new_path()
    post_df.to_csv(post_path, sep='\t', index=False)

    return pre_path, post_path
