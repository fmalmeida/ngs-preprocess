process PORECHOP {
  publishDir "${params.output}/preprocessing_outputs/nanopore/porechop", mode: 'copy'
  tag "${id}"
  label 'process_medium'

  input:
  file reads
  
  output:
  tuple val(id), path("${id}.trimmed.fq.gz"), val('nanopore')             optional true
  tuple val(id), path("${id}_porechop_barcodes/*.fq.gz"), val('nanopore') optional true

  when:
  !(reads =~ /input.*/)

  script:
  id = (reads.getBaseName() - ".fastq.gz" - ".fastq" - ".fq.gz" - ".fq")
  if (params.nanopore_is_barcoded)
  """
  # run porechop
  porechop -i ${reads} -t ${task.cpus} -b ${id}_porechop_barcodes --barcode_threshold 85

  # fix barcode extensions and gzip outputs
  cd ${id}_porechop_barcodes && \\
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
      -o ${id}.trimmed.fq

  # gzip output
  gzip ${id}.trimmed.fq
  """
}
