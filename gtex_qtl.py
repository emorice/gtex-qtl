"""
GTEx QTL calling pipeline
"""

import os

import wdl_galp

def broad_wdl(path):
    """
    Build path to the Broad's wdl scripts

    Args:
        path: path to the wanted wdl file relative to the root of the Broad
            repository
    """
    self_dir = os.path.dirname(os.path.realpath(__file__))

    return os.path.join(self_dir, 'gtex-pipeline', path)

# run the fastqtl tool through WDL
fastqtl = wdl_galp.run(broad_wdl('qtl/fastqtl.wdl'))

default_target = fastqtl
