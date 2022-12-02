"""
Preprocessing steps and utils.

By pre-processing, it is meant everything upstream of PEER or equivalent
correction
"""

try:
    import local_settings
except ModuleNotFoundError:
    pass

import os
import subprocess
import gzip
import contextlib

import galp
import wdl_galp
import pandas as pd

from .utils import broad_wdl

pbl = galp.Block()

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

    cov_df = extract_covariates(
            subject_path=subject_phenotypes_path,
            sample_path=sample_attributes_path
        )

    cov_df.to_csv(
            dst_path,
            sep='\t',
            index=False
        )
    return dst_path, list(cov_df.columns[1:])

@pbl.step
def gct_filter_columns(src_path, genotyped_subject_ids, _galp):
    """
    Drop the extra "id" column present in the public expression files, and
    subset samples with available genotype information.
    """
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

def extract_covariates(subject_path, sample_path):
    """
    Extract sex, WGS sequencing protocol and platform

    Returns:
        dataframe with columns SUBJID, SEX,
    """
    # Get subject id and sex
    subject_df = pd.read_table(subject_path, usecols=['SUBJID', 'SEX']).set_index('SUBJID')
    subject_df['SEX'] -= 1

    # get sample id and protocol
    sample_df = (
            # Keep sample ID, platform, and what sample was used for
            pd.read_table(sample_path, usecols=['SAMPID', 'SMGEBTCHT', 'SMAFRZE'])
            # Look only at samples used for WGS
            .query('SMAFRZE == "WGS"')
            .drop('SMAFRZE', axis=1)
            )

    protocol_coding = pd.DataFrame.from_records([
        {
            'SMGEBTCHT': 'PCR+ 30x Coverage WGS v2 (HiSeqX)',
            'WGS_PCR': 1,
            'WGS_HISEQX': 1,
            },
        {
            'SMGEBTCHT': 'PCR-Free 30x Coverage WGS v1 (HiSeqX)',
            'WGS_PCR': 0,
            'WGS_HISEQX': 1,
            },
        {
            'SMGEBTCHT': 'PCR+ 30x Coverage WGS (HiSeq2000)',
            'WGS_PCR': 1,
            'WGS_HISEQX': 0,
            },
        ], index='SMGEBTCHT')

    # define coded covariates for each sample
    sample_df = (
            sample_df
            .join(protocol_coding, on='SMGEBTCHT')
            .drop('SMGEBTCHT', axis=1)
            )

    # convert sample id to subject id
    sample_df.insert(0, 'SUBJID',
            'GTEX-' + sample_df['SAMPID'].str.split('-', expand=True)[1]
        )

    sample_df = ( sample_df
            .drop('SAMPID', axis=1)
            )

    # join and format
    cov_df = (
            sample_df
            .join(subject_df, on='SUBJID')
            .set_index('SUBJID')
            .T
            .reset_index()
            .rename({'index': 'ID'}, axis=1)
            )

    return cov_df

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
    Generate a lookup file matching samples to participants from a list of
    samples by simply interpreting the first two dash-separated components as
    the participant id. So "AA-BB-CC-DD" is considered to be a sample from
    participant "AA-BB".
    """
    path = _galp.new_path()
    with open(path, 'w', encoding='ascii') as stream:
        print('sample_id\tparticipant_id', file=stream)
        for sample_id in expression_sample_ids:
            participant_id = '-'.join(sample_id.split('-')[:2])
            print(f'{sample_id}\t{participant_id}', file=stream)
    return path

@pbl.step(vtag='0.2 extension')
def tabix(src_path, _galp):
    """
    Compress and index file with tabix
    """
    dst_path = _galp.new_path()

    # Try to preserve extension
    dst_path += '.'.join([''] + os.path.basename(src_path).split('.')[1:])
    dst_path += '.gz'

    with open(dst_path, 'wb') as stream:
        subprocess.run(['bgzip', '-c', src_path], stdout=stream, check=True)

    subprocess.run(['tabix', dst_path], check=True)
    index_path = dst_path + '.tbi'

    if not os.path.exists(index_path):
        raise RuntimeError('Tabix index not created')

    return dst_path, index_path

@pbl.step(vtag='0.2: both paths')
def indexed_vcf(_galp):
    """
    Symlink then index with tabix the vcf given in local_settings

    TODO: delete and replace with `tabix` step
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
def filtered_indexed_vcf(_galp, maf=0.01):
    """
    Filter out variants with less than `maf` maf, then recompress and index the
    result
    """
    dst_path = _galp.new_path() + '.vcf.gz'

    in_path = local_settings.GTEX_GENOTYPE_FILE

    with contextlib.ExitStack() as stack:
        vcf_p = stack.enter_context(subprocess.Popen(
            ['vcftools', '--gzvcf', in_path, '--maf', str(maf),
            '--stdout', '--recode', '--recode-INFO-all'],
            stdout=subprocess.PIPE))
        stream = stack.enter_context(open(dst_path, 'wb'))
        bgz_p = stack.enter_context(subprocess.Popen(['bgzip'],
            stdin=vcf_p.stdout, stdout=stream))
        vcf_p.stdout.close()
        if vcf_p.wait() or bgz_p.wait():
            raise RuntimeError('Vcf filtering failed')

    subprocess.run(['tabix', dst_path], check=True)

    return dst_path, dst_path + '.tbi'

@pbl.step
def chr_list(filtered_indexed_vcf, _galp):
    """
    File with the chromosome list from indexed vcf
    """
    path = _galp.new_path()
    with open(path, 'wb') as fobj:
        subprocess.run(['tabix', '-l', filtered_indexed_vcf[0]], stdout=fobj, check=True)
    return path

PREFIX='wb' # not really used but must be set to some string

def prepare_expression(wb_tpm, wb_counts, gene_model):
    """
    Pipeline meta-step to run the expression preparation step of the reference pipeline.

    All inputs are publicly available files. Code to retrieve them is provided
    in :mod:`gtex_qtl.downloads`.

    Args:
        wb_tpm: path to the raw expression (in TPMs)
        wb_counts: path to the raw expression (in read counts)
        gene_model: path to the gene model file of the reference pipeline.

    Returns:
        tuple: **expression_file** (path to the prepared expression) and
            **expression_file_index** (path to the matching index file)
    """
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

    return (
        prepared_expression[
                    'eqtl_prepare_expression_workflow'
                    '.eqtl_prepare_expression'
                    '.expression_bed'
                    ],
        prepared_expression[
                    'eqtl_prepare_expression_workflow'
                    '.eqtl_prepare_expression'
                    '.expression_bed_index'
                    ]
        )

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
