process filter {
  publishDir "${params.outdir}/longreads/${id}/filtered_reads", mode: 'copy'
  tag "filtering longreads with nanofilt"

  input:
  tuple val(id), file(reads), val(type)
  
  output:
  file "${custom_id}*"

  when:
  !(reads =~ /input.*/)

  script:
  if (params.nanopore_is_barcoded && type == 'nanopore') {
    custom_id = reads.getBaseName() - ".fastq.gz" - ".fastq"
  } else {
    custom_id = id
  }
  quality = (params.lreads_min_quality) ? "-q ${params.lreads_min_quality}" : ''
  length  = (params.lreads_min_length)  ? "-l ${params.lreads_min_length}" : ''
  """
  # Filtering
  gunzip -f -c $reads | NanoFilt ${quality} ${length} | gzip > ${custom_id}_filtered.fastq.gz ;
  """
}
