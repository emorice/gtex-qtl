"""
Tests for gtex_qtl.qtl
"""

import pandas as pd
import pytest

from gtex_qtl.qtl import call_qtls

VCF_CONTENT = """
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2
1\t10\tSNP1\tA\tG\t.\t.\t.\tGT\t0|0\t0|1
"""

@pytest.fixture
def vcf_path(tmp_path):
    """
    VCF fixture
    """
    path = tmp_path / 'test.vcf'
    with open(path, 'w', encoding='utf8') as fd:
        fd.write(VCF_CONTENT)
    return path

def test_call_qtls(vcf_path):
    """
    Returns well-formed calls
    """
    expression = pd.DataFrame({
        'chrom': ['1', '2', '2', '2', '3'],
        'start': [10, 10, 20, 30, 10],
        'sample1': [0.] * 5,
        })

    calls = call_qtls(
            (expression, 2),
            (20, 30),
            vcf_path,
            pd.DataFrame(), # gt covariates
            None # gx covariates
            )

    assert len(calls) == 3
    assert all(isinstance(item, pd.DataFrame) for item in calls)
