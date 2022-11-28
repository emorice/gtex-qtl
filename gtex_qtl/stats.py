"""
Generic statistical routines used throughout the pipelines
"""

import logging

import numpy as np
import scipy.optimize
from scipy.special import digamma, betaln

logger = logging.getLogger(__name__)

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

def fit_beta(mlp_mean, mlpc_mean, tol=1e-6, n_iter=1000):
    """
    Fit a beta distribution from sufficient statistics by maximum likelihood

    Uses very simple order-1 multiplicative updates, works usually well but not
    particularly fast.

    Args:
        mlp_mean: observed average negative log of values
            (:math:`\\frac 1 N \\sum{-\\ln(x_i)}`)
        mlpc_mean: observed average negative log of complementary values (one
            minus values, :math:`\\frac 1 N \\sum{-\\ln(1 - x_i)}`)
    """
    shape1, shape2 = 1., 1.

    i = 0
    mlp, mlpc = 1., 1.

    while (
        np.abs(mlp - mlp_mean) > tol * mlp_mean
        and np.abs(mlpc - mlpc_mean) > tol * mlpc_mean
        and i < n_iter
        ):
        shape1 *= mlp / mlp_mean
        shape2 *= mlpc / mlpc_mean
        dig1, dig2, dig12 = (
                digamma(shape1), digamma(shape2),
                digamma(shape1 + shape2)
                )
        mlp, mlpc = dig12 - dig1, dig12 - dig2
        i += 1

    if i >= n_iter:
        logger.warning(
            'Beta distribution fit did not converge in %d iterations', n_iter
            )
    else:
        logger.info(
            'Beta distribution fit converged in %d iterations', i
            )
    return shape1, shape2

def _log_scaled_t(values, dofs):
    return (
        - betaln(.5, .5 * dofs)
        - .5 * (dofs + 1) * np.log(1 + values**2)
        )

def fit_scaled_t(values):
    """
    Fit a scaled t distribution
    """
    opt = scipy.optimize.minimize_scalar(
            lambda param: - _log_scaled_t(values, np.exp(param)).sum()
            )
    if not opt.success:
        logger.warning(
                'Scaled-T distribution fit failed after %d iterations: %s',
                opt.nit, opt.message
            )
    else:
        logger.info(
                'Scaled-T distribution fit succeeded in %d iterations', opt.nit
            )
    return np.exp(opt.x)
