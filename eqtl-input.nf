#!/usr/bin/env nextflow

SAMPLES_GLOB = params.samples_glob
MAPPABILITY = params.mappability
METADATA = params.metadata
TPM = params.tpm
COUNTS = params.counts
GENOTYPE_PCA = params.genotype_pca
GTF = params.gtf
MAX_NUMBER_PCS = ['Lung': 300, 'Whole_blood': 400, 'PBMC': 300]

nextflow.enable.dsl=2

process make_tensorqtl_in {

	publishDir "${params.results}/tensorqtl-in", pattern: "*.tensorqtl-in.*"
    publishDir "${params.results}/pca", pattern: "*.PC-*"
    tag "${tissue}"
	memory '50 GB'
	time '168h'
	container 'docker://porchard/general:20230411143108'

	input:
	tuple val(tissue), path(samples), path(mappability), path(tpm), path(counts), path(metadata), path(gtf), path(genotype_pca)

	output:
	path("*.tsv")
	path("*.bed.gz")

	script:
	MAX_PCS = MAX_NUMBER_PCS.keySet().contains(tissue) ? MAX_NUMBER_PCS[tissue] : 100

	"""
	make-cis-eqtl-input.py --min-mappability 0.5 --mappability-scores $mappability --genotype-pcs 10 --min-phenotype-pcs 0 --max-phenotype-pcs $MAX_PCS --phenotype-pcs-step 5 $samples $tpm $counts $metadata $gtf $genotype_pca ${tissue}.
	"""

}


workflow {
	samples_for_scan = Channel.fromPath(SAMPLES_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, samples
    mappability = Channel.fromPath(MAPPABILITY)
    tpm = Channel.fromPath(TPM)
    counts = Channel.fromPath(COUNTS)
    metadata = Channel.fromPath(METADATA)
    gtf = Channel.fromPath(GTF)
    genotype_pca = Channel.fromPath(GENOTYPE_PCA)

    samples_for_scan.combine(mappability).combine(tpm).combine(counts).combine(metadata).combine(gtf).combine(genotype_pca) | make_tensorqtl_in
}
