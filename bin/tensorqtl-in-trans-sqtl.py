#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from qtl import norm, io
import argparse
import logging
from plutils import general, stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--genotype-pcs', default=10, dest='genotype_pcs', type=int)
parser.add_argument('--min-phenotype-pcs', default=0, dest='min_phenotype_pcs', type=int)
parser.add_argument('--max-phenotype-pcs', default=200, dest='max_phenotype_pcs', type=int)
parser.add_argument('--phenotype-pcs-step', default=5, dest='phenotype_pcs_step', type=int)
parser.add_argument('samples', help='List of sample IDs to include')
parser.add_argument('phenotypes', help='Intron excision matrix, using NWD IDs and already normalized')
parser.add_argument('metadata', help='')
parser.add_argument('genotype_pca', help='')
parser.add_argument('prefix')
args = parser.parse_args()

GENOTYPE_PCS = args.genotype_pcs
PHENOTYPE_PCS = list(range(args.min_phenotype_pcs, args.max_phenotype_pcs+1, args.phenotype_pcs_step))
SAMPLES = args.samples
PHENOTYPES = args.phenotypes
METADATA = args.metadata
GENOTYPE_PCA = args.genotype_pca
PREFIX = args.prefix

# SAMPLES = '/net/topmed10/working/porchard/rnaseq/work/eqtl-scan-in/Lung/samples-for-scan.txt'
# COUNT_MATRIX = '/net/topmed10/working/porchard/rnaseq/data/wolverine/all-cohorts.gene_counts.hdf5'
# GENOTYPE_PCS = 10
# GENOTYPE_PCA = '/net/topmed10/working/porchard/rnaseq/work/tensorqtl-in/lung/data/genotype-pcs.txt'
# PHENOTYPE_PCS = list(range(0, 20, 5))
# METADATA = '/net/topmed10/working/porchard/rnaseq/data/wolverine/metadata.txt'

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

logging.info('Loading sample list')
samples = pd.read_csv(SAMPLES, header=None)[0].astype(str).to_list()

logging.info('Loading metadata')
metadata = pd.read_csv(METADATA, sep='\t')
metadata.tor = metadata.tor.astype(str)
metadata = metadata[metadata.tor.isin(samples)]
tor_to_wgs = dict(zip(metadata.tor, metadata.wgs))

logging.info('Loading phenotypes')
phenotypes = pd.read_csv(PHENOTYPES, index_col=None, sep='\t')
SAMPLE_ORDER = phenotypes.columns.to_list()[4:]
assert(all(['NWD' in i for i in SAMPLE_ORDER]))
assert(pd.Series(SAMPLE_ORDER).value_counts().max() == 1)

logging.info('Inverse normalizing phenotypes')
invnorm_phenotypes = phenotypes.set_index(phenotypes.columns.to_list()[:4]).apply(norm.inverse_normal_transform, axis=1).reset_index()
assert(all(phenotypes.iloc[:,:4] == invnorm_phenotypes.iloc[:,:4]))
assert(invnorm_phenotypes.columns.to_list() == phenotypes.columns.to_list())
invnorm_phenotypes.to_csv(f'{PREFIX}tensorqtl-in.phenotypes.bed.gz', sep='\t', index=False)

# pca
logging.info('Performing PCA')
transformed, explained_variance, loadings = stats.pca(invnorm_phenotypes[SAMPLE_ORDER].T, return_loadings=True)
logging.info('Outputting PCA results')
loadings.to_csv(f'{PREFIX}PC-loadings.txt', sep='\t', index=True)
transformed.to_csv(f'{PREFIX}PC-scores.txt', sep='\t', index=True)
explained_variance = pd.DataFrame([['PC'+str(count), e] for count, e in enumerate(explained_variance, 1)], columns=['PC', 'fraction_variance_explained'])
explained_variance.to_csv(f'{PREFIX}PC-variance-explained.txt', sep='\t', index=False)

genotype_pcs = pd.read_csv(GENOTYPE_PCA, sep='\t', index_col=0).rename_axis(index='wgs').rename(columns=lambda x: f'genotype_{x}')
genotype_pcs = genotype_pcs.loc[SAMPLE_ORDER,[f'genotype_PC{i}' for i in range(1, GENOTYPE_PCS+1)]]

covariates = metadata[['inferred_sex', 'wgs', 'cohort']].rename(columns={'inferred_sex': 'sex'}).set_index('wgs')
covariates = covariates.loc[SAMPLE_ORDER,:]

if covariates.cohort.nunique() == 1:
    covariates = covariates.drop(columns=['cohort'])
covariates = pd.get_dummies(covariates, drop_first=True)
covariates = covariates.join(genotype_pcs)

phenotype_pca = transformed.rename(columns=lambda x: 'phenotype_' + x)
for phenotype_pcs in PHENOTYPE_PCS:
    if phenotype_pcs > len(phenotype_pca.columns):
        logging.warning(f"Skipping covariate file with {phenotype_pcs} phenotype PCs (there aren't that many)")
        continue
    covariates_with_phenotype_pcs = covariates.copy() if phenotype_pcs == 0 else covariates.join(phenotype_pca[[f'phenotype_PC{i}' for i in range(1, phenotype_pcs+1)]])
    covariates_with_phenotype_pcs.loc[SAMPLE_ORDER,:].T.to_csv(f'{PREFIX}tensorqtl-in.{phenotype_pcs}.covariates.tsv', sep='\t')    