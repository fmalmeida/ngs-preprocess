process NANOPACK {
  publishDir "${params.output}/preprocessing_outputs/${type}/QC", mode: 'copy'
  tag "${id}"
  label 'process_medium'

  input:
  tuple val(id), file(reads), val(type)
  
  output:
  path "${custom_id}*"

  when:
  !(reads =~ /input.*/)

  script:
  if (params.nanopore_is_barcoded && type == 'nanopore') {
    custom_id = reads.getBaseName() - ".fastq.gz" - ".fastq" - ".fq.gz" - ".fq"
  } else {
    custom_id = id
  }
  """
  # Plotting
  NanoPlot \\
      -t ${task.cpus} \\
      --fastq ${reads} \\
      -o ${custom_id}_nanoplot \\
      --N50 \\
      --title "${custom_id} sample" \\
      --plots hex dot kde ;

  # Checking Quality
  nanoQC \\
      -o ${custom_id}_nanoQC \\
      ${reads} ;

  # Generate Statistics Summary
  NanoStat \\
      --fastq ${reads} \\
      -t ${task.cpus} \\
      -n ${custom_id}.txt \\
      --outdir ${custom_id}_stats ;
  """
}
