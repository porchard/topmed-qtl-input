# RNA-seq paper scan input

## Dependencies

* Singularity (v. >=3)
* NextFlow (v. >= 21.04.0)


## Setup
 
1. Update nextflow.config as necessary for your compute environment.
2. Prep data: `make data`
3. Place the list of unrelated WGS samples here: `data/unrelated/unrelated.txt` (no header, a single column of IDs)
4. Place the genotype PCA here: `data/genotype-pca/PC-scores.txt` (first column: WGS ID, then columns PC1, PC2, ... PCN. Includes header)
5. Place metadata file here: `data/metadata/metadata.txt`
6. Place a file of pass-QC sample IDs here: `data/pass-qc/pass-qc.txt` (no header, a single column of RNA-seq sample IDs)
7. Place HDF5 files of gene expression here: `data/gene-expression/{counts,tpm}.hdf5`
8. [Run BAM remapping](https://github.com/porchard/topmed-rnaseq-remapping-and-leafcutter) and link repo here: `data/remapped` (ln -s /path/to/topmed-rnaseq-remapping-and-leafcutter data/remapped)

## Running

1. Select samples for cross-ancestry scans: `make eqtl-scan-in-lung eqtl-scan-in-PBMC eqtl-scan-in-monocyte eqtl-scan-in-T_cell eqtl-scan-in-whole-blood eqtl-scan-in-nasal-epithelial`
2. Make cis-eQTL and trans-eQTL input for cross-ancestry scans: `make tensorqtl-in-eQTL`
3. Generate leafcutter phenotypes: `make cluster-introns-no-nunique-mask`
4. Generate cis-sQTL and trans-sQTL scan input for cross-ancestry scans: `make tensorqtl-in-sQTL`