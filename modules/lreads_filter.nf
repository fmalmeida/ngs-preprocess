process FILTER {
  publishDir "${params.output}/final_output/${type}", mode: 'copy'
  tag "${id}"

  input:
  tuple val(id), file(reads), val(type)
  
  output:
  file "${custom_id}*"

  when:
  !(reads =~ /input.*/)

  script:
  if (params.nanopore_is_barcoded && type == 'nanopore') {
    custom_id = reads.getBaseName() - ".fastq.gz" - ".fastq" - ".fq.gz" - ".fq"
  } else {
    custom_id = id
  }
  quality = (params.lreads_min_quality) ? "-q ${params.lreads_min_quality}" : ''
  length  = (params.lreads_min_length)  ? "-l ${params.lreads_min_length}" : ''

  if (params.lreads_min_length || params.lreads_min_quality)
  """
  # filtering
  gunzip -f -c $reads | NanoFilt ${quality} ${length} | gzip > ${custom_id}.filtered.fq.gz ;
  """

  else
  """
  # save information that reads are not filtered
  cp $reads ${custom_id}.unfiltered.fq.gz
  """
}
