#!/usr/bin/env nextflow

nextflow.enable.dsl=2

COLLAPSED_GTF = params.collapsed_gtf
SAMPLES_GLOB = params.samples_glob // {tissue}..*, listing TOR IDs to used for each tissue
JUNCTION_FILES_GLOB = params.junction_files_glob
TOR_TO_NWD = params.tor_to_nwd

LEAFCUTTER_DIR = params.leafcutter_dir

process gather_junction_files {

    executor 'local'

    input:
    path("junction_files/*")

    output:
    path("junction_files")

    """
    echo pass
    """
}


process get_exon_list {

    container 'docker.io/porchard/general:20220406125608'
    memory '20 GB'
    cache 'lenient'

    input:
    path(collapsed_gtf)

    output:
    path('exons.txt.gz')

    """
    #!/usr/bin/env python
    import pandas as pd
    import qtl.annotation

    annot = qtl.annotation.Annotation('$collapsed_gtf')
    exon_df = pd.DataFrame([[g.chr, e.start_pos, e.end_pos, g.strand, g.id, g.name]
                            for g in annot.genes for e in g.transcripts[0].exons],
                        columns=['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name'])
    exon_df.to_csv('exons.txt.gz', sep='\t', index=False)
    """

}


process make_ratios {
    
    publishDir "${params.results}/ratios"
    container 'docker.io/porchard/leafcutter'
    memory '70 GB'
    cache 'lenient'
    time '168h'

    input:
    tuple val(tissue), path(samples), val(sample_size), path(junction_files_dir), path(exons), path(gtf), path(sample_participant_map)

    output:
    tuple val(tissue), path("${tissue}.leafcutter.bed.gz"), path("${tissue}.leafcutter.bed.gz.tbi"), path("${tissue}.leafcutter.phenotype_groups.txt")
    tuple val(tissue), path("${tissue}_perind.counts.gz"), emit: count_matrix
    path("${tissue}.leafcutter.PCs.txt")
    path("${tissue}_perind.counts.filtered.gz")
    path("${tissue}_perind_numers.counts.gz")
    path("${tissue}_pooled.gz")
    path("${tissue}_refined.gz")

    script:
    prefix = tissue
    min_clu_reads = Math.round(sample_size / 5)
    min_samples = Math.round(sample_size / 10)

    """
    ls ${junction_files_dir}/* | grep -f $samples > jf.txt
    cluster_prepare_fastqtl_use_custom_clustering_no_nunique_mask.py --leafcutter_dir ${LEAFCUTTER_DIR} --min_samples $min_samples --min_clu_reads $min_clu_reads jf.txt ${exons} ${gtf} ${prefix} ${sample_participant_map}
    """
}


process to_parquet {

    publishDir "${params.results}/parquet"
    memory { 30.GB * task.attempt }
    cache 'lenient'
    maxRetries 3
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}

    input:
    tuple val(tissue), path(x)

    output:
    tuple val(tissue), path("${tissue}.perind_counts.parquet")

    """
    sqtl-perind-counts-to-parquet.py $x ${tissue}.perind_counts.parquet
    """

}


process compute_masks {

    publishDir "${params.results}/masks"
    memory '100 GB'
    cache 'lenient'

    input:
    tuple val(tissue), path(x)

    output:
    tuple val(tissue), path("${tissue}.mask.parquet")

    """
    compute-sqtl-masks.py $x ${tissue}.mask.parquet
    """

}


process map_clusters_to_genes {

    publishDir "${params.results}/clusters-to-genes"
    memory '50 GB'
    cache 'lenient'
    container 'docker.io/porchard/leafcutter'

    input:
    tuple val(tissue), path(x), path(exons)

    output:
    tuple val(tissue), path("${tissue}.leafcutter.clusters_to_genes.txt")

    """
    Rscript /src/map_clusters_to_genes.R $x $exons ${tissue}.leafcutter.clusters_to_genes.txt
    """

}


process examine_samples_supporting_each_intron {

    publishDir "${params.results}/samples-supporting-each-intron"
    memory '50 GB'
    cache 'lenient'
    maxRetries 4
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}

    input:
    tuple val(tissue), path(samples), path(junction_files_dir)

    output:
    path("${tissue}.per-intron-support.txt")

    """
    ls ${junction_files_dir}/* | grep -f $samples > jf.txt
    examine-samples-supporting-each-intron.py jf.txt > ${tissue}.per-intron-support.txt
    """

}


workflow {
    gtf = Channel.fromPath(COLLAPSED_GTF)
    samples = Channel.fromPath(SAMPLES_GLOB).map({it -> [it.getName().tokenize('.')[0], it]})
    samples_with_sample_size = samples.map({it -> [it[0], it[1], it[1].countLines()]}) // tissue, file, N
    junctions = Channel.fromPath(JUNCTION_FILES_GLOB).toSortedList() | gather_junction_files
    tor_to_nwd = Channel.fromPath(TOR_TO_NWD)

    exon_list = get_exon_list(gtf)
    examine_samples_supporting_each_intron(samples.combine(junctions))
    
    ratios = samples_with_sample_size.combine(junctions).combine(exon_list).combine(gtf).combine(tor_to_nwd) | make_ratios
    to_parquet(ratios.count_matrix) | compute_masks
    map_clusters_to_genes(ratios.count_matrix.combine(exon_list))
}
