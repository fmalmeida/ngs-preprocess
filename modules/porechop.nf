process PORECHOP {
    publishDir "${params.output}/preprocessing_outputs/nanopore/porechop", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("**.fq.gz"), emit: fastqs

    when:
    meta.longreads_type == 'nanopore'

    script:
    if (params.nanopore_is_barcoded)
    """
    # run porechop
    porechop -i ${reads} -t ${task.cpus} -b ${meta.id}_porechop_barcodes --barcode_threshold 85

    # fix barcode extensions and gzip outputs
    cd ${meta.id}_porechop_barcodes && \\
        for i in *.fastq ; do mv \$i \${i%%.fastq}.fq ; done && \\
        for i in *.fq ; do gzip $i ; done
    """
    else
    """
    # run porechop
    porechop \\
        -i ${reads} \\
        -t ${task.cpus} \\
        --format fastq \\
        -o ${meta.id}.trimmed.fq

    # gzip output
    gzip ${meta.id}.trimmed.fq
    """
}
