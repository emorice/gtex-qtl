"""
Pipeline fork demonstrating pre-residualization of expression
"""

import numpy as np
import pandas as pd

import galp

pbl = galp.StepSet()

@pbl.step(vtag='0.1: no gzip')
def residualize(prepared_expression, combined_covariates, _galp):
    """
    Perform the residualization as a separate step and write back the result as
    a new expression file.

    Here residualization is taken as meaning regressing the prepared expression
    linearily on the covariates and taking the residuals of this regression.
    """

    # Shorthands: g - genes, s - samples, c - covariates

    # Read expression
    expr_df_gs = pd.read_table(prepared_expression)

    # Read covariates
    cov_df_cs = pd.read_table(combined_covariates)

    # Format as matching arrays
    samples = [ col for col in expr_df_gs if col.startswith('GTEX-') ]
    expr_gs = expr_df_gs[samples].to_numpy()
    cov_cs = cov_df_cs[samples].to_numpy()

    # Center covariates
    cov_cs -= cov_cs.mean(axis=1, keepdims=True)

    # Compute
    res_gs = (expr_gs - (
        (expr_gs @ cov_cs.T) @ np.linalg.inv(cov_cs @ cov_cs.T) @ cov_cs
        ))

    # Reformat as dataframe
    res_df_gs = pd.DataFrame(res_gs, columns=samples)

    # Reinclude other columns
    res_df_gs = pd.concat((
        expr_df_gs.drop(samples, axis=1),
        res_df_gs
        ), axis=1)

    # Write back as uncompressed tsv
    dst_path = _galp.new_path() + '.bed'
    res_df_gs.to_csv(dst_path, sep='\t', index=False)

    return dst_path
