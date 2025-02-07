#!/usr/bin/env python
# coding: utf-8


import sys
import pandas as pd
import numpy as np

PARQUET_COUNT_MATRIX = sys.argv[1]
OUT = sys.argv[2]

mat = pd.read_parquet(PARQUET_COUNT_MATRIX)


def get_masks(counts_df):
    calculate_frac = lambda x: float(x[0])/float(x[1]) if x[1] > 0 else 0
    frac_df = counts_df.applymap(lambda x: calculate_frac([int(i) for i in x.split('/')]))
    pct_zero = (frac_df == 0).sum(1) / frac_df.shape[1]  # for zero counts, frac is zero
    n_unique = frac_df.apply(lambda x: len(x.unique()), axis=1)
    zscore_df = ((frac_df.T-frac_df.mean(1)) / frac_df.std(1)).T # per-individual z-scores w/in each feature

    n = np.floor(frac_df.shape[1]*0.1)
    if n < 10:
        n = 10
    pass_pct_zero_mask = (pct_zero <= 0.5)
    pass_n_unique_mask = (n_unique >= n)
    # additional filter for low complexity
    ns = zscore_df.shape[1]
    zscore_mask = ((zscore_df.abs()<0.25).sum(1) >= ns-3) & ((zscore_df.abs() > 6).sum(1) <= 3) # first term: nearly all z-scores have small magnitude (little variability); second term: very few z-scores have large magnitude
    pass_zscore_mask = ~zscore_mask

    assert(all(pass_pct_zero_mask.index == pass_n_unique_mask.index))
    assert(all(pass_pct_zero_mask.index == pass_zscore_mask.index))

    mask_df = pd.DataFrame({'pass_pct_zero_mask': pass_pct_zero_mask, 'pass_n_unique_mask': pass_n_unique_mask, 'pass_zscore_mask':pass_zscore_mask })
    return mask_df


masked = get_masks(mat)
masked.to_parquet(OUT)
