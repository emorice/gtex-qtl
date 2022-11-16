"""
Reference PEER correction
"""

import wdl_galp

from .utils import broad_wdl

def peer_and_combine(expression_file, gpcs_covariates, additional_covariates):
    """
    Reference pipeline meta-step, compute peer factors and combine all covariates

    Args:
        expression_file: path to file with prepared expression
        gpcs_covariates: path to file with genotype principal components
        additional_covariates: path to file with other covariates to include
    Returns:
        path to combined covariates file
    """
    combined_covariates = wdl_galp.run(broad_wdl('qtl/eqtl_peer_factors.wdl'),
            ** { f'eqtl_peer_factors.{key}': value
                for key, value in {
                    'expression_file': expression_file,
                    'prefix': 'tissue', # placeholder, we don't use the prefixes
                    'num_peer': 60,  # "60 factors for N ≥ 350"
                    # optional inputs
                    'genotype_pcs': gpcs_covariates,
                    'add_covariates': additional_covariates,
                    # runtime
                    'memory': 0,
                    'disk_space': 0,
                    'num_threads': 64,
                    'num_preempt': 0,
                    }.items()
                }
            )

    return combined_covariates[
                'eqtl_peer_factors_workflow'
                '.eqtl_peer_factors'
                '.combined_covariates'
                ]
