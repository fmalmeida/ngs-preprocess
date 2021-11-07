process nanopack {
  publishDir "${params.outdir}/longreads/nanopack_out", mode: 'copy'
  tag "Checking longreads qualities with Nanopack"
  //validExitStatus 0,1 // To momentainaly fix problem with matplotlib

  input:
    file reads
  output:
    file "${id}*"

  script:
  id = (reads.getBaseName() - "fastq.gz" - ".fastq")
  """
  # Plotting
  NanoPlot -t ${params.threads} --fastq ${reads} -o ${id}_nanoplot --N50 --title "${id} sample" --plots hex dot pauvre kde ;

  # Checking Quality
  nanoQC -o ${id}_nanoQC ${reads} ;

  # Generate Statistics Summary
  NanoStat --fastq ${reads} -t ${params.threads} -n ${id}.txt --outdir ${id}_stats ;
  """
}
