"""
GTEx QTL calling pipeline

This file defines the reference reproduction pipeline. Sister modules contain
pipeline variants.
"""

import os
import re
import gzip
import contextlib
import urllib.request
import shutil
import subprocess

import pandas as pd
import plotly.graph_objects as go

import galp
import wdl_galp

from . import utils
from . import residualize
from . import compare

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

def local_wdl(path):
    """
    Build path to our modified or new wdl scripts
    """
    # Package dir itself
    self_dir = os.path.dirname(os.path.realpath(__file__))

    return os.path.join(self_dir, path)


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

_GTEX_BASE_URL = 'https://storage.googleapis.com/gtex_analysis_v8/'

wb_tpm = gct_filter_columns(urlretrieve(_GTEX_BASE_URL +
        'rna_seq_data/gene_tpm/gene_tpm_2017-06-05_v8_whole_blood.gct.gz'))

wb_counts = gct_filter_columns(urlretrieve(_GTEX_BASE_URL +
        'rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_whole_blood.gct.gz'))

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
expression_file_index = prepared_expression[
            'eqtl_prepare_expression_workflow'
            '.eqtl_prepare_expression'
            '.expression_bed_index'
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

combined_covariates_file = combined_covariates[
            'eqtl_peer_factors_workflow'
            '.eqtl_peer_factors'
            '.combined_covariates'
            ]

# 4) Run fastqtl
# ==============

def run_fastqtl(expression_files, covariates_file=None):
    """
    Meta step to build the fastqtl task
    """
    return wdl_galp.run(local_wdl('fastqtl.wdl'),
        expression_bed=expression_files[0],
        expression_bed_index=expression_files[1],
        vcf=indexed_vcf[0], vcf_index=indexed_vcf[1],
        prefix=PREFIX,
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

fastqtl = run_fastqtl(
        (expression_file, expression_file_index),
        covariates_file=combined_covariates_file
        )

# 5) Get published results and compare
# ====================================


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
def reference_tissue_expression(reference_expression, ref_tissue_name):
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

expression_shapes_cmp = {
        'reference': expression_shape(reference_tissue_expression),
        'computed': expression_shape(expression_file)
        }

reference_covariates = untar(urlretrieve(_GTEX_BASE_URL +
        'single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_covariates.tar.gz'
        ))

pbl.bind(reference_results=untar(urlretrieve(_GTEX_BASE_URL +
    'single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar'
    )))

@pbl.step
def reference_tissue_egenes(reference_results, ref_tissue_name):
    """
    Path to specific tissue results
    """
    return os.path.join(reference_results,
        'GTEx_Analysis_v8_eQTL',
        f'{ref_tissue_name}.v8.egenes.txt.gz'
        )

pbl.bind(computed_tissue_egenes=fastqtl['fastqtl_workflow.fastqtl_permutations_merge.genes'])

@pbl.view
def qvalues_cmp(reference_tissue_egenes, computed_tissue_egenes):
    """
    Scatter plot of log-q-values
    """
    merged = (
        pd.read_table(reference_tissue_egenes)
        .merge(
            pd.read_table(computed_tissue_egenes),
            on='gene_id',
            suffixes=['_reference', '_computed']
            )
        )

    return go.Figure(
        data=[
            go.Scatter(
                x=merged['qval_reference'],
                y=merged['qval_computed'],
                mode='markers',
                marker={'size': 3},
                name='recomputed'
                ),
            go.Scatter(
                x=merged['qval_reference'],
                y=merged['qval_reference'],
                mode='lines',
                name='y = x'
                )
            ],
        layout={
            'title': 'Comparison of published and reproduced gene statistics',
            'xaxis': {'type': 'linear', 'title': 'Gene-level q-value (published)'},
            'yaxis': {
                'type': 'linear',
                'title': 'Gene-level q-value (recomputed)'
                },
            'height': 1000,
            }
        )

# 6) Pre-residualization alternatives
# ===================================

residualized_expression = tabix(residualize.residualize(expression_file,
        combined_covariates_file))

residualized_fastqtl = run_fastqtl(residualized_expression)

pbl.bind(res_vs_orig_raster=compare.datashader_scatter(
        compare.all_pvals(fastqtl),
        compare.all_pvals(residualized_fastqtl),
        log=True
        ))

@pbl.view
def residualized_pvals_plot(res_vs_orig_raster):
    """
    Residualized vs original p-values for all pairs
    """
    fig = compare.plot_ds_scatter(res_vs_orig_raster)
    fig.update_layout({
        'title':
        'Comparison of nominal p-values for all tested gene-variant pairs',
        'xaxis.title': 'Reproduced p-value',
        'yaxis.title': 'Pre-residualized p-value',
        })
    return fig

blind_linear_expression = tabix(
        residualize.residualize_blind_linear(
            expression_file,
            gpcs_covariates,
            additional_covariates[0],
            )
        )

blind_linear_fastqtl = run_fastqtl(blind_linear_expression)

pbl.bind(blind_vs_res_raster=compare.datashader_scatter(
        compare.all_pvals(residualized_fastqtl),
        compare.all_pvals(blind_linear_fastqtl),
        log=True
        ))

@pbl.view
def blind_pvals_plot(blind_vs_res_raster):
    """
    Blind linear vs residualized p-values for all pairs
    """
    fig = compare.plot_ds_scatter(blind_vs_res_raster)
    fig.update_layout({
        'title':
        'Comparison of nominal p-values for all tested gene-variant pairs',
        'xaxis.title': 'Pre-residualized p-value',
        'yaxis.title': 'Blind linear p-value',
        })
    return fig

# END
# ===

plots = [ residualized_pvals_plot, blind_pvals_plot ]

default_target = plots
