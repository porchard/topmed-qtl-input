#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging
import argparse

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

#METADATA = '/lab/work/porchard/daily/2021-06-14/metadata.txt'
#UNRELATED = '/lab/work/porchard/daily/2021-06-14/kingunrelated.txt'
#MATCHES = '/lab/work/porchard/daily/2021-06-14/best-matches.txt'
#METADATA, UNRELATED, MATCHES, PASS_PCA_QC, TISSUE = sys.argv[1:]

parser = argparse.ArgumentParser()
parser.add_argument('--metadata')
parser.add_argument('--cohorts', nargs='+')
parser.add_argument('--tissues', nargs='+')
parser.add_argument('--unrelated')
parser.add_argument('--pass-pca-qc', dest='pass_pca_qc', help='List of TOR IDs that passed QC')
parser.add_argument('--in-genotype-pca', dest='in_genotype_pca', help='List of samples in the genotype PCA')
args = parser.parse_args()

METADATA = args.metadata
TISSUES = args.tissues
COHORTS = args.cohorts
UNRELATED = args.unrelated
PASS_PCA_QC = args.pass_pca_qc
IN_GENOTYPE_PCA = args.in_genotype_pca


SAMPLE_DROPOUT = [] # cohort, stage, stage_name, n

unrelated = pd.read_csv(UNRELATED, sep='\t', header=None)[0].to_list()
pass_pca_qc = pd.read_csv(PASS_PCA_QC, sep='\t', header=None)[0].to_list()
in_genotype_pca = pd.read_csv(IN_GENOTYPE_PCA, sep='\t', header=None)[0].to_list()

metadata = pd.read_csv(METADATA, sep='\t', index_col=False)
metadata = metadata[metadata.tissue.isin(TISSUES)]
metadata = metadata[metadata.cohort.isin(COHORTS)]

metadata = metadata[metadata.has_rnaseq==1]
for cohort in COHORTS:
    SAMPLE_DROPOUT.append([cohort, 1, 'with_rnaseq', (metadata.cohort==cohort).sum()])

metadata = metadata[metadata.tor.isin(pass_pca_qc)]
for cohort in COHORTS:
    SAMPLE_DROPOUT.append([cohort, 2, 'not_pca_outlier', (metadata.cohort==cohort).sum()])

metadata = metadata[~metadata.wgs.isnull()]
for cohort in COHORTS:
    SAMPLE_DROPOUT.append([cohort, 3, 'has_wgs_match', (metadata.cohort==cohort).sum()])

metadata = metadata[metadata.wgs.isin(unrelated)]
for cohort in COHORTS:
    SAMPLE_DROPOUT.append([cohort, 4, 'in_unrelated_set', (metadata.cohort==cohort).sum()])

metadata = metadata[metadata.inferred_sex.isin(['female', 'male'])]
for cohort in COHORTS:
    SAMPLE_DROPOUT.append([cohort, 5, 'successfully_inferred_sex', (metadata.cohort==cohort).sum()])

# if the same WGS sample appears multiple times, randomly select one (except favor certain timepoints for MESA/SPIROMICS)
def select_rnaseq_sample_from_subject(df):
    tmp = df.copy()
    if tmp.cohort.nunique() != 1:
        logging.warning('Multiple cohorts: {} ({})'.format(tmp.wgs.unique()[0], ', '.join(tmp.cohort.unique())))
    if tmp.cohort.unique()[0] == 'SPIROMICS':
        if 'BASELINE' in tmp.collection_visit.to_list():
            tmp = tmp[tmp.collection_visit=='BASELINE']
    if tmp.cohort.unique()[0] == 'MESA':
        if 'MESA_exam5' in tmp.collection_visit.to_list():
            tmp = tmp[tmp.collection_visit=='MESA_exam5']
        elif 'MESA_exam1' in tmp.collection_visit.to_list():
            tmp = tmp[tmp.collection_visit=='MESA_exam1']
    return tmp.sample(frac=1.0, random_state=2021).tor.to_list()[0]

metadata = metadata.sample(frac=1.0, random_state=2021).reset_index(drop=True)
USE_FOR_SCAN = [select_rnaseq_sample_from_subject(df) for subject, df in metadata.groupby('wgs')]
metadata = metadata[metadata.tor.isin(USE_FOR_SCAN)]
for cohort in COHORTS:
    SAMPLE_DROPOUT.append([cohort, 6, 'after_selecting_one_rnaseq_per_wgs', (metadata.cohort==cohort).sum()])

if not all(metadata.wgs.isin(in_genotype_pca)):
    metadata = metadata[metadata.wgs.isin(in_genotype_pca)]
    for cohort in COHORTS:
        SAMPLE_DROPOUT.append([cohort, 7, 'in_genotype_pca', (metadata.cohort==cohort).sum()])

# plot sample dropout
dropout = pd.DataFrame(SAMPLE_DROPOUT, columns=['cohort', 'stage', 'stage_name', 'n'])
dropout.stage_name = pd.Categorical(dropout.stage_name, categories=dropout.stage_name.unique(), ordered=True)
dropout.to_csv('sample-dropout.tsv', sep='\t', index=False)

labels = {k: '{} (n={:,})'.format(k, v) for k, v in dropout.groupby('stage_name').n.sum().to_dict().items()}
dropout['stage_name'] = dropout.stage_name.map(labels)

cohort_labels = {k: '{}\n(final n={:,})'.format(k, v) for k, v in dropout.groupby('cohort').n.min().items()}
dropout.cohort = dropout.cohort.map(cohort_labels)

fig, ax = plt.subplots(figsize=(10,6))
sns.barplot(x='cohort', y='n', hue='stage_name', data=dropout)
ax.set_xlabel('')
ax.set_ylabel('RNA-seq samples')
ax.legend().set_title('')
fig.tight_layout()
fig.savefig('sample-dropout.png', dpi=300)
fig.clf()

# output
metadata.sort_values('cohort').loc[:,['tor']].to_csv(sys.stdout, sep='\t', index=False, header=False)