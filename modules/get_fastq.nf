process GET_FASTQ {
  publishDir "${params.output}/SRA_FETCH/FASTQ", mode: 'copy'
  tag "${sra_ids.replaceAll(~/\s/,'')}"
  label 'process_low'
  
  input:
  val(sra_ids)

  output:
  path "*_data", emit: sra_fastq

  when: sra_ids

  script:
  def sra_ids = "${sra_ids.replaceAll(~/\s/,'')}"
  """
  fasterq-dump \\
    --include-technical \\
    --split-files \\
    --threads $task.cpus \\
    --outdir ./${sra_ids}_data \\
    --progress \\
    $sra_ids
  """
}
