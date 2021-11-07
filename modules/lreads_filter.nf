process filter {
  publishDir "${params.outdir}/longreads/filtered_reads", mode: 'copy'
  tag "filtering longreads with nanofilt"

  input:
    file reads
  output:
    file "${id}*"

  script:
  id = (reads.getBaseName() - "fastq.gz" - ".fastq")
  quality = (params.lreads_min_quality) ? "-q ${params.lreads_min_quality}" : ''
  length  = (params.lreads_min_length)  ? "-l ${params.lreads_min_length}" : ''
  """
  # Filtering
  gunzip -f -c $reads | NanoFilt ${quality} ${length} | gzip > ${id}_filtered.fastq.gz ;
  """
}
