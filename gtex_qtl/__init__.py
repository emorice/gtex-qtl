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

import pandas as pd

import galp
import wdl_galp

from . import utils

try:
    import local_settings
except ModuleNotFoundError:
    pass

# pylint: disable=redefined-outer-name

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

# 0.1) Extract additional covariates and genotyped subject list
# =============================================================

pbl.bind(subject_phenotypes_path=local_settings.GTEX_SUBJECT_PHENOTYPES_FILE)
pbl.bind(sample_attributes_path=local_settings.GTEX_SAMPLE_ATTRIBUTES_FILE)

@pbl.step(vtag='0.5: naming', items=2)
def additional_covariates(subject_phenotypes_path, sample_attributes_path, _galp):
    """
    Generate a covariate file with extra covariates.

    Also returns the list of the available subject IDs. Availability mostly
    means that a WGS sample was found in the metadata files for a given subject.
    """
    dst_path = _galp.new_path()

    cov_df = utils.extract_covariates(
            subject_path=subject_phenotypes_path,
            sample_path=sample_attributes_path
        )

    cov_df.to_csv(
            dst_path,
            sep='\t',
            index=False
        )
    return dst_path, list(cov_df.columns[1:])

# 0.2) Download public expression files and gene model
# ====================================================

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

@pbl.step
def gct_filter_columns(src_path, additional_covariates, _galp):
    """
    Drop the extra "id" column present in the public expression files, and
    subset samples with available genotype information.
    """
    genotyped_subject_ids = additional_covariates[1]

    dst_path = _galp.new_path() + '.gct.gz'
    with gzip.open(src_path, 'rt', encoding='ascii') as src:
        with gzip.open(dst_path, 'wt', encoding='ascii') as dst:
            # Headers
            dst.write(src.readline())
            dst.write(src.readline())

            # Table
            expr_df = pd.read_table(src)

            sample_ids = [
                sample_id
                for sample_id in expr_df.columns
                if sample_id.startswith('GTEX-')
                if f"GTEX-{sample_id.split('-')[1]}" in genotyped_subject_ids
                ]

            expr_df[["Name", "Description", *sample_ids]].to_csv(
                    dst, sep='\t', index=False
                    )
    return dst_path

_GTEX_BASE_URL = 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/'

wb_tpm = gct_filter_columns(urlretrieve(_GTEX_BASE_URL +
        'gene_tpm/gene_tpm_2017-06-05_v8_whole_blood.gct.gz'))

wb_counts = gct_filter_columns(urlretrieve(_GTEX_BASE_URL +
        'gene_reads/gene_reads_2017-06-05_v8_whole_blood.gct.gz'))

_GTEX_GENE_MODEL_URL = (
    'https://personal.broadinstitute.org/francois/topmed/'
    'gencode.v26.GRCh38.ERCC.genes.gtf.gz'
    )

gene_model = urlretrieve(_GTEX_GENE_MODEL_URL, gunzip=True)


# 0.3) Prepare lookup table and chromosome list
# =============================================

@pbl.step
def expression_sample_ids(counts):
    """
    List sample IDs from expression file
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
def sample_participant_lookup(expression_sample_ids, _galp):
    """
    Generate a lookup file for all sample IDs found in expression file
    """
    path = _galp.new_path()
    with open(path, 'w', encoding='ascii') as stream:
        print('sample_id\tparticipant_id', file=stream)
        for sample_id in expression_sample_ids:
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
    with open(path, 'wb') as fobj:
        subprocess.run(['tabix', '-l', indexed_vcf[0]], stdout=fobj, check=True)
    return path

# 1) Prepare expression (filter and normalize)
# ============================================

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
                # Runtime section (not really used)
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

# 1.1) Prepare genotype PCs
# =========================

@pbl.step
def filter_format_pcs(genotype_pcs_path, _galp):
    """
    Filter first PCs and format them into a covariate file
    """
    pcs_df = pd.read_table(genotype_pcs_path,
            usecols = ['FID'] + [ f'PC{i+1}' for i in range(5) ]
            )

    pcs_df['SUBJID'] = 'GTEX-' + pcs_df['FID'].str.split('-', expand=True)[1]
    pcs_df = (pcs_df
        .drop('FID', axis=1)
        .set_index('SUBJID')
        .T
        .reset_index()
        .rename({'index': 'ID'}, axis=1)
        )

    dst_path = _galp.new_path()
    pcs_df.to_csv(dst_path, sep='\t', index=False)
    return dst_path

gpcs_covariates = filter_format_pcs(local_settings.GTEX_GENOTYPE_PCS_FILE)

# 2) PEER factors and 3) combine covariates
# =========================================

combined_covariates = wdl_galp.run(broad_wdl('qtl/eqtl_peer_factors.wdl'),
        ** { f'eqtl_peer_factors.{key}': value
            for key, value in {
                'expression_file': expression_file,
                'prefix': PREFIX,
                'num_peer': 60,  # "60 factors for N ≥ 350"
                # optional inputs
                'genotype_pcs': gpcs_covariates,
                'add_covariates': additional_covariates[0],
                # runtime
                'memory': 0,
                'disk_space': 0,
                'num_threads': 64,
                'num_preempt': 0,
                }.items()
            }
        )

# 4) Run fastqtl
# ==============

fastqtl = wdl_galp.run(broad_wdl('qtl/fastqtl.wdl'),
        expression_bed=expression_file,
        expression_bed_index=prepared_expression[
            'eqtl_prepare_expression_workflow'
            '.eqtl_prepare_expression'
            '.expression_bed_index'],
        vcf=indexed_vcf[0], vcf_index=indexed_vcf[1],
        prefix=PREFIX,
        )

default_target = combined_covariates # fastqtl
