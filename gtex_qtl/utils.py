"""
Covariate extraction util
"""

import pandas as pd

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
