"""
GTEx QTL calling pipeline
"""

import os
import urllib.request
import shutil

import galp
import wdl_galp

pbl = galp.StepSet()

def broad_wdl(path):
    """
    Build path to the Broad's wdl scripts

    Args:
        path: path to the wanted wdl file relative to the root of the Broad
            repository
    """
    self_dir = os.path.dirname(os.path.realpath(__file__))

    return os.path.join(self_dir, 'gtex-pipeline', path)

# download public expression files
@pbl.step(vtag='0.2: return path')
def urlretrieve(url, _galp):
    """
    Download file from url and check it in the store
    """
    to_path = _galp.new_path()
    with urllib.request.urlopen(url) as in_fd:
        with open(to_path, 'wb') as out_fd:
            shutil.copyfileobj(in_fd, out_fd)
    return to_path

_GTEX_BASE_URL = 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/'

wb_tpm = urlretrieve(_GTEX_BASE_URL +
        'gene_tpm/gene_tpm_2017-06-05_v8_whole_blood.gct.gz')

wb_counts = urlretrieve(_GTEX_BASE_URL +
        'gene_reads/gene_reads_2017-06-05_v8_whole_blood.gct.gz')

_GTEX_GENE_MODEL_URL = (
    'https://personal.broadinstitute.org/francois/topmed/'
    'gencode.v26.GRCh38.ERCC.genes.gtf.gz'
    )

gene_model = urlretrieve(_GTEX_GENE_MODEL_URL)

# various normalization steps
prepared_expression = wdl_galp.run(broad_wdl('qtl/eqtl_prepare_expression.wdl'),
        **{
            'eqtl_prepare_expression.tpm_gct': wb_tpm,
            'eqtl_prepare_expression.counts_gct': wb_counts,
            'eqtl_prepare_expression.annotation_gtf': gene_model,
            }
        )

# run the fastqtl tool through WDL
fastqtl = wdl_galp.run(broad_wdl('qtl/fastqtl.wdl'),
        expression_bed=prepared_expression['expression_bed']
        )

default_target = fastqtl
