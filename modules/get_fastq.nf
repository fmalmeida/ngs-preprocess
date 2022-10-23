process GET_FASTQ {
  publishDir "${params.output}/SRA_FETCH/PREFETCH", mode: 'copy'
  tag "$sra_ids"
  label 'process_medium'
  
  input:
  val(sra_ids)

  output:
  path "${sra_ids}_data", emit: sra_fastq

  script:
  """
  fastq-dump --gzip --split-3 --outdir ./${sra_ids}_data $sra_ids
  """
}
