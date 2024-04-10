#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import gzip
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

#import glob
#JUNCTION_FILES = glob.glob('/net/topmed3/working/porchard/rnaseq/work/intron-cluster-exploration-subset-whole-blood/work/0b/855a44a497c5cc7a5affbf7ea9d71b/junction_files/TOR100*')
JUNCTION_FILES_LIST = sys.argv[1]
MAX_INTRON_LENGTH = 500000
CHROMS = set([f'chr{i}' for i in range(1, 23)] + ['chrX'])

JUNCTION_FILES = pd.read_csv(JUNCTION_FILES_LIST, header=None)[0].to_list()

counts = dict() # intron --> [n_samples, total_counts]
for file_count, f in enumerate(JUNCTION_FILES, 1):
    sample = os.path.basename(f).split('.')[0]
    logging.info('Processing file {} ({:,} of {:,})'.format(f, file_count, len(JUNCTION_FILES)))
    with gzip.open(f, 'rt') as fh:
        for line in fh:
            chrom, start, end, name, score, strand, rA, rb, rgb, blockCount, blockSize, blockStarts = line.rstrip().split()
            if chrom not in CHROMS:
                continue
            if strand not in ['+', '-']:
                continue
            Aoff, Boff = blockSize.split(",")
            intron_start, intron_end = (int(start) + int(Aoff), int(end) - int(Boff) + 1)
            score = int(score)
            if intron_end - intron_start > MAX_INTRON_LENGTH:
                continue
            intron = f'{chrom}:{intron_start}:{intron_end}:{strand}'
            if intron not in counts:
                counts[intron] = [0, 0]
            counts[intron][0] += 1
            counts[intron][1] += score


counts = [[k] + v for k, v in counts.items()]
df = pd.DataFrame(counts, columns=['intron', 'n_samples', 'n_reads'])
df.to_csv(sys.stdout, sep='\t', index=False)