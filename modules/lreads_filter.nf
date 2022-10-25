process FILTER {
    publishDir "${params.output}/final_output/${meta.longreads_type}", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    
    output:
    path "*filtered.fq.gz", emit: fastqs

    when:
    !(reads =~ /input.*/)

    script:
    quality = (params.lreads_min_quality) ? "-q ${params.lreads_min_quality}" : ''
    length  = (params.lreads_min_length)  ? "-l ${params.lreads_min_length}"  : ''

    if (params.lreads_min_length || params.lreads_min_quality)
    """
    # filtering
    gunzip -f -c $reads | NanoFilt ${quality} ${length} | gzip > ${meta.id}.filtered.fq.gz ;
    """

    else
    """
    # save information that reads are not filtered
    cp $reads ${meta.id}.unfiltered.fq.gz
    """
}
