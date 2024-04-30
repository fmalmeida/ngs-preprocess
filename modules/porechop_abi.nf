process PORECHOP_ABI {
    publishDir "${params.output}/preprocessing_outputs/nanopore/porechop_abi", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("**.fq.gz"), emit: fastqs

    when:
    meta.longreads_type == 'nanopore'

    script:
    """
    # run porechop
    porechop_abi \\
        --ab_initio \\
        -i ${reads} \\
        -o ${meta.id}.trimmed.fq

    # gzip output
    gzip ${meta.id}.trimmed.fq
    """
}
