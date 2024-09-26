"""
Tests for gtex_qtl.qtl
"""

import subprocess
import numpy as np
import pandas as pd
import pytest

from gtex_qtl.qtl import call_qtls

@pytest.fixture
def vcf(tmp_path):
    """
    VCF fixture
    """
    path = tmp_path / 'test.vcf'
    rng = np.random.default_rng(0)
    n_samples = 30
    n_sites = 10
    gts = rng.choice(2, size=(n_sites, n_samples, 2))
    with open(path, 'w', encoding='utf8') as fd:
        ## Header
        # std columns
        fd.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
        # samples names
        for i in range(n_samples):
            fd.write(f'\tsample{i+1}')
        fd.write('\n')

        ## Data
        for i_site, gt_site in enumerate(gts):
            # genotype metadata
            fd.write(f'1\t{i_site+1}\tSNP{i_site}\tA\tG\t.\t.\t.\tGT')
            # genotype calls
            for gt in gt_site:
                fd.write(f'\t{gt[0]}|{gt[1]}')
            fd.write('\n')

    subprocess.check_call(f'bgzip {path}', shell=True)
    return f'{path}.gz', gts

def test_call_qtls(vcf) -> None:
    """
    Returns well-formed calls
    """
    vcf_path, gts = vcf
    rng = np.random.default_rng(1)
    n_samples = 30
    n_genes = 1000

    # Draw a common random factor
    factor_s = rng.normal(size=n_samples)
    # Draw genes with large share of variance from factor
    expr_sg = 0.97 * factor_s[:, None] + 0.22 * rng.normal(size=(n_samples, n_genes))
    #  Add a large genetic effect to first ten genes
    expr_sg[:, :10] += .5 * gts.sum(-1).T

    expression = pd.DataFrame({
        # One gene on test chrom, all others as controls on other chrom
        'chrom': ['1'] * (n_genes // 2) + ['2'] * (n_genes // 2),
        'start': np.arange(n_genes) + 1,
        } | {
            f'sample{i+1}': expr_g
            for i, expr_g in enumerate(expr_sg)
        }
        )

    covariates = pd.DataFrame({
        'variable': ['sex'],
        } | {
            f'sample{i+1}': [rng.choice(2)]
            for i in range(n_samples)
        })

    gt_regressors = [
            { 'method': 'external', 'data': covariates},
            { 'method': 'auto', 'data': None}
            ]
    gx_regressors = [
            { 'method': 'external', 'data': covariates},
            { 'method': 'auto', 'data': None}
            ]

    calls = call_qtls(
            (expression, 2),
            (0, 10),
            vcf_path,
            qtl_config={
                'num_null_genes': 450,
                'genotype_regressors': gt_regressors, # gt covariates
                'expression_regressors': gx_regressors, # gx covariates
                },
            )

    assert len(calls) == 3
    assert all(isinstance(item, pd.DataFrame) for item in calls)
    for item in calls:
        print(item.to_string())
