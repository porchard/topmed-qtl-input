#!/usr/bin/env python
# coding: utf-8


import sys
import pandas as pd

INPUT, OUTPUT = sys.argv[1:]


ne = pd.read_csv(INPUT, sep=' ', dtype='str')
ne = ne.set_index('chrom')

ne = ne.rename(columns=lambda x: x.replace('.regtools_junc.txt.gz', ''))
ne.to_parquet(OUTPUT)