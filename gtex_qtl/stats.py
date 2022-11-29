"""
Generic statistical routines used throughout the pipelines
"""

import logging

import numpy as np
import scipy.optimize
import scipy.special
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

def _log_convergence(opt, task):
    if not opt.success:
        logger.warning(
                '%s failed after %d iterations: %s', task,
                opt.nit, opt.message
            )
    else:
        logger.info(
                '%s succeeded in %d iterations', task, opt.nit
            )

def fit_scaled_t(values):
    """
    Fit a scaled t distribution
    """
    opt = scipy.optimize.minimize_scalar(
            lambda param: - _log_scaled_t(values, np.exp(param)).sum()
            )
    _log_convergence(opt, 'Scaled-T distribution fit')
    return np.exp(opt.x)

def st_log_pdf(st2, dofs):
    """
    Log pdf of a scaled-t variable.

    A Scaled-T variable is one distributed like X / sqrt(nu) when X is a
    student-T with nu degrees of freedom.

    Args:
        st2: square of the observed scaled-t
        dofs: degrees of freedom
    """
    return (
        - .5 * (dofs + 1) * np.log(1 + st2)
        - scipy.special.betaln(0.5, 0.5 * dofs)
    )

def st_sf(st2, dofs):
    """
    Survival function of a scaled-t variable. See :func:`st_log_pdf`
    """
    return scipy.special.betainc(.5 * dofs, .5, 1. / (1. + st2))


def mle_log_rep_max(sfs):
    """
    Given observations of the maxima of iid, identically-sized random vectors,
    estimate size of the original vectors.

    If X_ik, 1 <= k <= K, are iid with survival function G, and we observe  M_i
    = max(X_ik, 1 <= k <= K), we can estimate K by maximum likelihood provided
    we know G(M_i). Intuitively, the larger K is, the further in the tail of G
    we are, and the smaller the values of G(M).

    Args:
        sfs: values of the survival function of the original random variable at
            each max.

    Returns:
        log of the maximum likelihood estimator of the number of observation
        maxed over
    """
    return (
        np.log(len(sfs))
        - np.log(-np.sum(scipy.special.log1p(-sfs)))
    )

def _maxed_scaled_t_ll(st2s, dofs):
    """
    Log-likelihood of the max of several scaled-t observations
    """
    return (
        np.sum(st_log_pdf(st2s, dofs))
        + len(st2s) * mle_log_rep_max(st_sf(st2s, dofs))
    )

def fit_max_scaled_t(best_scaled_t2s):
    """
    Fit max-of-many scaled-t distribution

    Args:
        best_scaled_t2s: observed maxima of the squared scaled ts.

    Returns:
        tuple with the mle of the number of degrees of freedom, and the number
        of observations maxed over.
    """
    opt = scipy.optimize.minimize_scalar(
            lambda dofs: -_maxed_scaled_t_ll(best_scaled_t2s, dofs),
            (0.1, 0.2))
    _log_convergence(opt, 'Max-Scaled-T distribution fit')
    mle_dofs = opt.x
    return (
            mle_dofs,
            np.exp(mle_log_rep_max(st_sf(best_scaled_t2s, mle_dofs)))
            )

def _fq_t_loss(st2s, dofs):
    """
    Re-implementation of the loss function used by FastQTL to re-estimate DoFs
    """
    pvalues = st_sf(st2s, dofs)
    mean = np.mean(pvalues)
    var = np.var(pvalues, ddof=1)
    concentration_est = mean * (1. - mean) / var - 1
    alpha_est = mean * concentration_est
    return np.abs(alpha_est - 1.)

def fqtl_fit_max_scaled_t(best_scaled_t2s):
    """
    Re-implementation of the method used by FastQTL to re-estimate DoFs and
    multiple p-value distribution
    """

    # Step 1: find the number of degrees of freedom that makes the
    # method-of-moment estimate of the p-value beta shape parameter 1 closest to
    # 1
    fq_opt = scipy.optimize.minimize_scalar(
            lambda dofs: _fq_t_loss(best_scaled_t2s, dofs),
            (0.1, 0.2))
    _log_convergence(fq_opt, 'FQTL estimation of DoFs')

    est_dofs = fq_opt.x

    # Step 2: fit a beta distribution by mle to pvalues computed with such DoFs
    pvalues = st_sf(best_scaled_t2s, est_dofs)
    beta_shapes = fit_beta(
            -np.mean(np.log(pvalues)),
            -np.mean(np.log(1. - pvalues))
            )

    return est_dofs, beta_shapes
