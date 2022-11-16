"""
Reference FastQTL analysis
"""

import wdl_galp

from .utils import local_wdl

def run_fastqtl(expression_files, indexed_vcf, gene_model, covariates_file=None):
    """
    Meta step to build the fastqtl task

    Args:
        expression_files: tuple with path to file and path to index
        indexed_vcf: tuple with path to vcf and path to index
        gene_model: path to gene model file
    """
    return wdl_galp.run(local_wdl('fastqtl.wdl'),
        expression_bed=expression_files[0],
        expression_bed_index=expression_files[1],
        vcf=indexed_vcf[0], vcf_index=indexed_vcf[1],
        prefix='tissue', # placeholder as we don't need prefixes
        permutations="1000 10000", # from readme and 2020 methods
        chunks=100, # from readme
        fdr=0.05, # from 2020 methods
        annotation_gtf=gene_model,
        # optional parameters
        maf_threshold=0.01, # from 2020 methods
        **( {'covariates': covariates_file}
            if covariates_file else {}),
        # runtime parameters
        **{
            f'{step}.{key}': value
            for key, value in {
                    'memory': 0,
                    'disk_space': 0,
                    'num_threads': 1,
                    'num_preempt': 0,
                }.items()
            for step in [
                'fastqtl_nominal',
                'fastqtl_permutations_scatter',
                'fastqtl_permutations_merge',
                'fastqtl_postprocess'
                ]
            }
        )
