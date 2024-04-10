#!/usr/bin/env nextflow

SAMPLES_GLOB = params.samples_glob
PHENOTYPES_GLOB = params.phenotypes_glob
METADATA = params.metadata
GENOTYPE_PCA = params.genotype_pca

nextflow.enable.dsl=2

process make_tensorqtl_in {

	publishDir "${params.results}/tensorqtl-in", pattern: "*.tensorqtl-in.*"
    publishDir "${params.results}/pca", pattern: "*.PC-*"
	container 'docker.io/porchard/general:20220406125608'
    tag "${tissue}"
	memory '50 GB'
	maxForks 5
	executor 'local'

	input:
	tuple val(tissue), path(samples), path(leafcutter_phenotypes), path(metadata), path(genotype_pca)

	output:
	path("*.tsv")
	path("${tissue}.tensorqtl-in.phenotypes.bed.gz")

	"""
	tensorqtl-in-trans-sqtl.py --genotype-pcs 15 --min-phenotype-pcs 0 --max-phenotype-pcs 50 --phenotype-pcs-step 1 $samples $leafcutter_phenotypes $metadata $genotype_pca ${tissue}.
	"""

}


workflow {
	samples_for_scan = Channel.fromPath(SAMPLES_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, samples
	phenotypes = Channel.fromPath(PHENOTYPES_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, phenotypes
    metadata = Channel.fromPath(METADATA)
    genotype_pca = Channel.fromPath(GENOTYPE_PCA)

    samples_for_scan.combine(phenotypes, by: 0).combine(metadata).combine(genotype_pca) | make_tensorqtl_in
}
