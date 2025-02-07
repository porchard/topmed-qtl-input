#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from qtl import norm, io
import argparse
import logging
from plutils import general, stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--mappability-scores', dest='mappability_scores', default=None, help='File of gene mappability scores (optional unless mappability threshold imposed). Mappability filter imposed only during PCA')
parser.add_argument('--min-mappability', default=None, dest='min_mappability', type=float)
parser.add_argument('--genotype-pcs', default=10, dest='genotype_pcs', type=int)
parser.add_argument('--min-phenotype-pcs', default=0, dest='min_phenotype_pcs', type=int)
parser.add_argument('--max-phenotype-pcs', default=200, dest='max_phenotype_pcs', type=int)
parser.add_argument('--phenotype-pcs-step', default=5, dest='phenotype_pcs_step', type=int)
parser.add_argument('samples', help='List of sample IDs to include')
parser.add_argument('tpm', help='TPM HDF5 matrix')
parser.add_argument('counts', help='counts HDF5 matrix')
parser.add_argument('metadata', help='')
parser.add_argument('gtf', help='GTF file')
parser.add_argument('genotype_pca', help='')
parser.add_argument('prefix')
args = parser.parse_args()

MAPPABILITY = args.mappability_scores
MIN_MAPPABILITY = args.min_mappability
GENOTYPE_PCS = args.genotype_pcs
PHENOTYPE_PCS = list(range(args.min_phenotype_pcs, args.max_phenotype_pcs+1, args.phenotype_pcs_step))
SAMPLES = args.samples
TPM_MATRIX = args.tpm
COUNT_MATRIX = args.counts
METADATA = args.metadata
GTF = args.gtf
GENOTYPE_PCA = args.genotype_pca
PREFIX = args.prefix

KEEP_CHROMS = [f'chr{i}' for i in list(range(1, 23)) + ['X']]

# SAMPLES = '/net/topmed10/working/porchard/rnaseq/work/eqtl-scan-in/Lung/samples-for-scan.txt'
# TPM_MATRIX = '/net/topmed10/working/porchard/rnaseq/data/wolverine/all-cohorts.gene_tpm.hdf5'
# COUNT_MATRIX = '/net/topmed10/working/porchard/rnaseq/data/wolverine/all-cohorts.gene_counts.hdf5'
# MAPPABILITY = '/net/topmed10/working/porchard/rnaseq/data/mappability/hg38_gene_mappability.txt.gz'
# MIN_MAPPABILITY = 0.5
# GTF = '/net/topmed10/working/porchard/rnaseq/work/tensorqtl-in/lung/flattened.gtf'
# GENOTYPE_PCS = 10
# GENOTYPE_PCA = '/net/topmed10/working/porchard/rnaseq/work/tensorqtl-in/lung/data/genotype-pcs.txt'
# PHENOTYPE_PCS = list(range(0, 20, 5))
# METADATA = '/net/topmed10/working/porchard/rnaseq/data/wolverine/metadata.txt'

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

logging.info('Loading gene locations')
phenotype_locations = io.gtf_to_tss_bed(GTF)
phenotype_locations.gene_id = phenotype_locations.gene_id
phenotype_locations = phenotype_locations[phenotype_locations.chr.isin(KEEP_CHROMS)]


mappability = None
if MAPPABILITY is not None:
    logging.info('Loading mappability scores')
    mappability = pd.read_csv(MAPPABILITY, header=None, sep='\t', names=['gene', 'mappability'])
    mappability['gene_no_version'] = mappability.gene.map(lambda x: general.strip_version_from_ensembl_gene_id(x, keep_par_y=True))
    mappability = dict(zip(mappability.gene_no_version, mappability.mappability))


logging.info('Loading sample list')
samples = pd.read_csv(SAMPLES, header=None)[0].astype(str).to_list()


logging.info('Loading TPM')
tpm = pd.read_hdf(TPM_MATRIX)
tpm = tpm.loc[samples,:]
logging.info('Generating list of genes that pass TPM threshold')
pass_tpm_filter_in_fraction_samples = (tpm>0.1).astype(int).mean()
genes_passing_tpm_filter = pass_tpm_filter_in_fraction_samples[pass_tpm_filter_in_fraction_samples>=0.2].index.to_list()
logging.info('{:,} of {:,} genes pass threshold'.format(len(genes_passing_tpm_filter), len(tpm.columns)))


logging.info('Loading counts')
counts = pd.read_hdf(COUNT_MATRIX)
counts = counts.loc[samples,:]

logging.info('Filtering to chromosomes: {}'.format(', '.join(KEEP_CHROMS)))
counts = counts[list(set(counts.columns).intersection(set(phenotype_locations.gene_id)))]

logging.info('Performing TMM normalization of counts')
tmm = norm.edger_cpm(counts.T).T


logging.info('Filtering to genes that passed TPM filter')
tmm_filtered = tmm[list(set(genes_passing_tpm_filter).intersection(set(phenotype_locations.gene_id)))].T

logging.info('Performing inverse normal transformation')
inv_norm = norm.inverse_normal_transform(tmm_filtered).T

logging.info('Outputting inverse normal transformed gene counts')
inv_norm.T.to_csv(f'{PREFIX}inverse-normalized-gene-expression.csv', index=True, sep='\t')




if mappability is not None:
    logging.info('Performing mappability filter')
    m = inv_norm.columns.to_series().map(lambda x: general.strip_version_from_ensembl_gene_id(x, keep_par_y=True)).map(mappability)
    DROP = m[(m.isnull()) | (m<MIN_MAPPABILITY)].index.to_list()
    logging.info('Keeping {:,} of {:,} genes'.format(len(m) - len(DROP), len(m)))
    for g in DROP:
        logging.info(f'Dropping gene {g} for PCA')
    inv_norm_mappability_filtered = inv_norm.loc[:,~inv_norm.columns.isin(DROP)]


# pca
logging.info('Performing PCA')
transformed, explained_variance, loadings = stats.pca(inv_norm_mappability_filtered, return_loadings=True)
logging.info('Outputting PCA results')
loadings.to_csv(f'{PREFIX}PC-loadings.txt', sep='\t', index=True)
transformed.to_csv(f'{PREFIX}PC-scores.txt', sep='\t', index=True)
explained_variance = pd.DataFrame([['PC'+str(count), e] for count, e in enumerate(explained_variance, 1)], columns=['PC', 'fraction_variance_explained'])
explained_variance.to_csv(f'{PREFIX}PC-variance-explained.txt', sep='\t', index=False)


logging.info('Loading metadata')
metadata = pd.read_csv(METADATA, sep='\t')
metadata.tor = metadata.tor.astype(str)
metadata = metadata[metadata.tor.isin(samples)]


tor_to_wgs = dict(zip(metadata.tor, metadata.wgs))
phenotypes = inv_norm.T.rename_axis(index='gene_id').rename(columns=tor_to_wgs)
assert(all(phenotypes.columns.to_series().str.contains('NWD')))
assert(phenotypes.columns.to_series().value_counts().max() == 1)

phenotypes = phenotype_locations.join(phenotypes, how='inner')


genotype_pcs = pd.read_csv(GENOTYPE_PCA, sep='\t', index_col=0).rename_axis(index='wgs').reset_index()
genotype_pcs = genotype_pcs[genotype_pcs.wgs.isin(phenotypes.columns)].set_index('wgs')
genotype_pcs.columns = ['genotype_{}'.format(i) for i in genotype_pcs.columns]
genotype_pcs = genotype_pcs[[f'genotype_PC{i}' for i in range(1, args.genotype_pcs+1)]]
# verify that all samples appear in the genotype PCA
assert(len([i for i in samples if tor_to_wgs[i] not in set(genotype_pcs.index)]) == 0)


covariates = metadata[['inferred_sex', 'wgs', 'cohort']].rename(columns={'inferred_sex': 'sex'})
covariates = covariates[covariates.wgs.isin(phenotypes.columns.to_list())]
covariates = covariates.set_index('wgs')
if covariates.cohort.nunique() == 1:
    covariates = covariates.drop(columns=['cohort'])
covariates = pd.get_dummies(covariates, drop_first=True)
covariates = covariates.join(genotype_pcs)

# now output the files
SAMPLE_ORDER = phenotypes.columns.to_list()[4:]
phenotypes.rename(columns={'chr': '#chr'}).to_csv(f'{PREFIX}tensorqtl-in.phenotypes.bed.gz', sep='\t', index=False)

phenotype_pca = transformed.rename(columns=lambda x: 'phenotype_' + x)
phenotype_pca.index = phenotype_pca.index.map(tor_to_wgs)
for phenotype_pcs in PHENOTYPE_PCS:
    if phenotype_pcs > len(phenotype_pca.columns):
        logging.warning(f"Skipping covariate file with {phenotype_pcs} phenotype PCs (there aren't that many)")
        continue
    covariates_with_phenotype_pcs = covariates.copy() if phenotype_pcs == 0 else covariates.join(phenotype_pca[[f'phenotype_PC{i}' for i in range(1, phenotype_pcs+1)]])
    covariates_with_phenotype_pcs.loc[SAMPLE_ORDER,:].T.to_csv(f'{PREFIX}tensorqtl-in.{phenotype_pcs}.covariates.tsv', sep='\t')
