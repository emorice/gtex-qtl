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
    subject_df = pd.read_table(subject_path, usecols=['SUBJID', 'SEX']).set_index('SUBJID')

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

    sample_df = (
            sample_df
            .join(protocol_coding, on='SMGEBTCHT')
            .drop('SMGEBTCHT', axis=1)
            )

    sample_df.insert(0, 'SUBJID',
            'GTEX-' + sample_df['SAMPID'].str.split('-', expand=True)[1]
        )

    sample_df = ( sample_df
            .drop('SAMPID', axis=1)
            )

    return sample_df.join(subject_df, on='SUBJID')
