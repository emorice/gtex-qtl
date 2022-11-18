"""
Generic statistical routines used throughout the pipelines
"""

import numpy as np

def regress(data_bs, covariates_cs):
    """
    Regress data linearly against covariates, adding an intercept.

    Args:
        data_bs: response variables to residualize as numby array with a
            (**b**) batch dimension then the (**s**) samples dimension matching the
            covariates
        covariates_cs: covariates to use as features, with the first dimensions
            going through said covariates and the second the (**s**) samples
            matching the response data.
    """
    # Add the intercept pseudo-covariate
    covariates_cs = np.vstack((np.ones(data_bs.shape[-1]), covariates_cs))

    # Compute
    residuals_bs = (data_bs - (
        (data_bs @ covariates_cs.T)
        @ np.linalg.inv(covariates_cs @ covariates_cs.T)
        @ covariates_cs
        ))

    return residuals_bs

def regress_missing(data_bs, valid_bs, covariates_cs):
    """
    Variant of :func:`regress` that takes a validity mask on data. Invalid cells
    of `data_bs` must be zeroed.

    This is **not** equivalent to regressing individually each vector of data
    after filtering samples, as this would require solving a linear system with
    a different set of samples for each one.
    """
    # Add the intercept pseudo-covariate
    covariates_cs = np.vstack((np.ones(data_bs.shape[-1]), covariates_cs))

    # Compute only one covariate variance covariance estimator, with all samples
    _n_vars, n_samples = covariates_cs.shape
    cov_cc = (covariates_cs @ covariates_cs.T) / n_samples

    # Compute the individual cov-gt covariances with variable samples
    cov_cb = (covariates_cs @ data_bs.T) / np.sum(valid_bs, -1)

    weights_cb = np.linalg.inv(cov_cc) @ cov_cb

    return data_bs - weights_cb.T @ covariates_cs
