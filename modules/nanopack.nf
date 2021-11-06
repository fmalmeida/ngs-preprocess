process nanopack {
  publishDir "${params.outdir}/longreads/nanopack_out", mode: 'copy'
  tag "Checking longreads qualities with Nanopack"
  //validExitStatus 0,1 // To momentainaly fix problem with matplotlib

  input:
    file reads
    val threads
  output:
    file "${id}*"

  script:
  id = (reads.getBaseName() - "fastq.gz" - ".fastq")
  """
  # Plotting
  #NanoPlot -t $threads --fastq ${reads} -o ${id}_nanoplot --N50 --title "${id} sample" --plots hex dot pauvre kde ;

  # nanoplot is currently having problems with pauvre ... let's not do it
  NanoPlot -t $threads --fastq ${reads} -o ${id}_nanoplot  --N50 --title "${id} sample" --plots hex dot kde ;

  # Checking Quality
  nanoQC -o ${id}_nanoQC ${reads} ;

  # Generate Statistics Summary
  NanoStat --fastq ${reads} -n ${id}.txt --outdir ${id}_stats ;
  """
}
