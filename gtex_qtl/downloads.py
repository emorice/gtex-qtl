"""
Public data retrieval
"""

import os
import re
import contextlib
import gzip
import urllib
import shutil

import galp
import pandas as pd

pbl = galp.Block()

@pbl.step(vtag='0.3: suffixes')
def urlretrieve(url, _galp, preserve_suffix=True, gunzip=False):
    """
    Download file from url and check it in the store
    """

    if preserve_suffix:
        # Try to preserve the file extension as a dowstream step may rely on it
        # for file typing
        suffix = '.' + '.'.join(url.split('/')[-1].split('.')[1:])
        if gunzip:
            suffix = re.sub('.gz$', '', suffix)
    else:
        suffix = ''
    to_path = _galp.new_path() + suffix

    with contextlib.ExitStack() as stack:
        in_fd = stack.enter_context(urllib.request.urlopen(url))
        if gunzip:
            in_fd = stack.enter_context(gzip.open(in_fd))
        out_fd = stack.enter_context(open(to_path, 'wb'))
        shutil.copyfileobj(in_fd, out_fd)
    return to_path

_GTEX_BASE_URL = 'https://storage.googleapis.com/gtex_analysis_v8/'
_GTEX_GENE_MODEL_URL = (
    'https://personal.broadinstitute.org/francois/topmed/'
    'gencode.v26.GRCh38.ERCC.genes.gtf.gz'
    )

input_files = dict(
    wb_tpm = urlretrieve(_GTEX_BASE_URL +
            'rna_seq_data/gene_tpm/gene_tpm_2017-06-05_v8_whole_blood.gct.gz'),
    wb_counts = urlretrieve(_GTEX_BASE_URL +
            'rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_whole_blood.gct.gz'),
    gene_model = urlretrieve(_GTEX_GENE_MODEL_URL, gunzip=True)
    )
"""
Publicly available files needed as inputs of the pipeline.

Note that the expression data files cannot be used as-is.
    * They contain an extra ID column that is unexpected by the rest of the pipeline
    * They contain samples for which no genotype is available.

In :mod:`gtex_qtl.preprocess` utils to filter out the extra columns are provided.
"""

@pbl.step
def untar(path, _galp):
    """
    Extract tar archive
    """
    dst_path = _galp.new_path()
    os.makedirs(dst_path, exist_ok=True)
    shutil.unpack_archive(path, dst_path)
    return dst_path

pbl.bind(reference_expression=untar(urlretrieve(_GTEX_BASE_URL +
        'single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_expression_matrices.tar'
        )))

pbl.bind(ref_tissue_name='Whole_Blood')

@pbl.step
def published_tissue_expression(reference_expression, ref_tissue_name):
    """
    Path to specific tissue in extracted archive

    Tissue name is capitalized, underscore separated
    """
    return os.path.join(reference_expression,
            'GTEx_Analysis_v8_eQTL_expression_matrices',
            f'{ref_tissue_name}.v8.normalized_expression.bed.gz'
            )

@pbl.step
def expression_shape(expression_file):
    """
    Extract basic counts from expression file
    """
    expr_df = pd.read_table(expression_file)

    return {
        'num_genes': len(expr_df),
        'num_samples': sum(col.startswith('GTEX-') for col in expr_df)
        }

pbl.bind(reference_results=untar(urlretrieve(_GTEX_BASE_URL +
        'single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar'
    )))

reference_covariates = untar(urlretrieve(_GTEX_BASE_URL +
        'single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_covariates.tar.gz'
        ))

@pbl.step
def published_tissue_egenes(reference_results, ref_tissue_name):
    """
    Path to specific tissue results
    """
    return os.path.join(reference_results,
        'GTEx_Analysis_v8_eQTL',
        f'{ref_tissue_name}.v8.egenes.txt.gz'
        )
