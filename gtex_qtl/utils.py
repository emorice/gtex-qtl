"""
File finders
"""

import os

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

def local_wdl(path):
    """
    Build path to our modified or new wdl scripts
    """
    # Package dir itself
    self_dir = os.path.dirname(os.path.realpath(__file__))

    return os.path.join(self_dir, path)
