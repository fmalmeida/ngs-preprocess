process nanopack {
  publishDir "${params.outdir}/longreads/${id}/nanopack_out", mode: 'copy'
  tag "Checking longreads qualities with Nanopack"
  //validExitStatus 0,1 // To momentainaly fix problem with matplotlib

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
  """
  # Plotting
  NanoPlot -t ${params.threads} --fastq ${reads} -o ${custom_id}_nanoplot --N50 --title "${custom_id} sample" --plots hex dot kde ;

  # Checking Quality
  nanoQC -o ${custom_id}_nanoQC ${reads} ;

  # Generate Statistics Summary
  NanoStat --fastq ${reads} -t ${params.threads} -n ${custom_id}.txt --outdir ${custom_id}_stats ;
  """
}
