"""
GTEx QTL calling pipeline
"""

import os
import re
import gzip
import contextlib
import urllib.request
import shutil
import subprocess

import galp
import wdl_galp

from . import utils

try:
    import local_settings
except ModuleNotFoundError:
    pass

pbl = galp.StepSet()

def broad_wdl(path):
    """
    Build path to the Broad's wdl scripts

    Args:
        path: path to the wanted wdl file relative to the root of the Broad
            repository
    """
    # Containing directory
    self_dir = os.path.dirname(
            # Package dir itself
            os.path.dirname(os.path.realpath(__file__))
                )

    return os.path.join(self_dir, 'gtex-pipeline', path)

# download public expression files
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

@pbl.step(vtag='0.2: file extension')
def gct_drop_id(src_path, _galp):
    """
    Drop the extra "id" column present in the public expression files
    """
    dst_path = _galp.new_path() + '.gct.gz'
    with gzip.open(src_path, 'rt', encoding='ascii') as src:
        with gzip.open(dst_path, 'wt', encoding='ascii') as dst:
            # Headers
            dst.write(src.readline())
            dst.write(src.readline())

            # Table head
            line = src.readline()
            if not line.startswith('id\t'):
                raise TypeError('No extra id column at beginning of table')
            dst.write(line.split('\t', maxsplit=1)[1])

            # Table body
            for line in src:
                dst.write(line.split('\t', maxsplit=1)[1])

    return dst_path

_GTEX_BASE_URL = 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/'

wb_tpm = gct_drop_id(urlretrieve(_GTEX_BASE_URL +
        'gene_tpm/gene_tpm_2017-06-05_v8_whole_blood.gct.gz'))

wb_counts = gct_drop_id(urlretrieve(_GTEX_BASE_URL +
        'gene_reads/gene_reads_2017-06-05_v8_whole_blood.gct.gz'))

_GTEX_GENE_MODEL_URL = (
    'https://personal.broadinstitute.org/francois/topmed/'
    'gencode.v26.GRCh38.ERCC.genes.gtf.gz'
    )

gene_model = urlretrieve(_GTEX_GENE_MODEL_URL, gunzip=True)

@pbl.step
def sample_ids(counts):
    """
    Extract sample IDs from expression file
    """
    with gzip.open(counts, 'rt') as stream:
        stream.readline() # Version tag
        stream.readline() # Table shape
        header = stream.readline()

        ids = header.strip().split()

        # Name, Description, [samples]
        generic_ids, gtex_ids = ids[:2], ids[2:]

        if any(gen_id.startswith('GTEX') for gen_id in generic_ids):
            raise TypeError('Expected 2 non-GTEX columns at start of file')
        if not all(gtex_id.startswith('GTEX') for gtex_id in gtex_ids):
            raise TypeError('Expected only GTEX ids after third column')

        return gtex_ids

@pbl.step
def sample_participant_lookup(sample_ids, _galp):
    """
    Generate a lookup file for participant IDs
    """
    path = _galp.new_path()
    with open(path, 'w') as stream:
        print('sample_id\tparticipant_id', file=stream)
        for sample_id in sample_ids:
            participant_id = '-'.join(sample_id.split('-')[:2])
            print(f'{sample_id}\t{participant_id}', file=stream)
    return path

@pbl.step(vtag='0.2: both paths')
def indexed_vcf(_galp):
    """
    Symlink then index with tabix the vcf given in local_settings
    """
    path = _galp.new_path()
    os.symlink(
        local_settings.GTEX_GENOTYPE_FILE,
        path
        )
    subprocess.run(['tabix', path], check=True)
    index_path = path + '.tbi'

    if not os.path.exists(index_path):
        raise RuntimeError('Tabix index not created')

    return path, index_path

@pbl.step
def chr_list(indexed_vcf, _galp):
    """
    File with the chromosome list from indexed vcf
    """
    path = _galp.new_path()
    with open(path, 'w') as fobj:
        subprocess.run(['tabix', '-l', indexed_vcf[0]], stdout=fobj, check=True)
    return path

# 1) various normalization steps
PREFIX='wb' # not really used but must be set to some string

prepared_expression = wdl_galp.run(broad_wdl('qtl/eqtl_prepare_expression.wdl'),
        **{ f'eqtl_prepare_expression.{key}': value
            for key, value in {
                'tpm_gct': wb_tpm,
                'counts_gct': wb_counts,
                'annotation_gtf': gene_model,
                'sample_participant_ids':
                    sample_participant_lookup(counts=wb_counts),
                'vcf_chr_list': chr_list,
                'prefix': PREFIX,
                # Runtime section
                'memory': 0,
                'disk_space': 0,
                'num_threads': 0,
                'num_preempt': 0,
            }.items()
        }
    )
expression_file = prepared_expression[
            'eqtl_prepare_expression_workflow'
            '.eqtl_prepare_expression'
            '.expression_bed'
            ]

# 2) PEER factors
peer_factors = wdl_galp.run(broad_wdl('qtl/eqtl_peer_factors.wdl'),
        ** { f'eqtl_peer_factors.{key}': value
            for key, value in {
                'expression_file': expression_file,
                'prefix': PREFIX,
                'num_peer': 60,  # "60 factors for N ≥ 350"
                # runtime
                'memory': 0,
                'disk_space': 0,
                'num_threads': 64,
                'num_preempt': 0,
                }.items()
            }
        )

# 3) Genotype PCs and combine covariates

@pbl.step(vtag='0.3: fix order')
def extract_covariates(subject_path, sample_path, _galp):
    """
    Generate a covariate file with extra covariates
    """
    dst_path = _galp.new_path()
    utils.extract_covariates(
            subject_path=subject_path,
            sample_path=sample_path
        ).to_csv(
            dst_path,
            sep='\t',
            index=False
        )
    return dst_path

additional_covariates = extract_covariates(
        subject_path=local_settings.GTEX_SUBJECT_PHENOTYPES_FILE,
        sample_path=local_settings.GTEX_SAMPLE_ATTRIBUTES_FILE
        )

# 4) run fastqtl
fastqtl = wdl_galp.run(broad_wdl('qtl/fastqtl.wdl'),
        expression_bed=expression_file,
        expression_bed_index=prepared_expression[
            'eqtl_prepare_expression_workflow'
            '.eqtl_prepare_expression'
            '.expression_bed_index'],
        vcf=indexed_vcf[0], vcf_index=indexed_vcf[1],
        prefix=PREFIX,
        )

default_target = peer_factors # fastqtl
