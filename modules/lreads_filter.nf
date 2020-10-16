process lreads_filter {
  publishDir "${params.outdir}/longreads/filtered_reads", mode: 'copy'
  // Loads the necessary Docker image
  container 'fmalmeida/ngs-preprocess'
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
  source activate nanopack;

  # Filtering
  gunzip -f -c $reads | NanoFilt ${quality} ${length} | gzip > ${id}-filtered-reads.fastq.gz ;
  """
}
