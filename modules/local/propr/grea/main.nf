process PROPR_GREA {
    tag "$meta.id"
    label 'process_high'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-limma_r-ggplot2_r-propr:209490acb0e524e3' :
        'community.wave.seqera.io/library/bioconductor-limma_r-ggplot2_r-propr:17abd3f137436739' }"

    input:
    tuple val(meta), path(adj)
    tuple val(meta2), path(gmt)

    output:
    tuple val(meta), path("*.grea.tsv"), emit: results
    path "versions.yml",                 emit: versions
    path "*.R_sessionInfo.log",          emit: session_info

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'grea.R'
}
