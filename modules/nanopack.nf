process nanopack {
  publishDir "${params.outdir}/longreads", mode: 'copy'
  // Loads the necessary Docker image
  container 'fmalmeida/ngs-preprocess'
  tag "Checking longreads qualities with Nanopack"

  input:
    file reads
    val threads
  output:
    file "${reads.baseName}*"

  script:
    """
  source activate nanopack;
  # Plotting
  NanoPlot -t $threads --fastq ${reads} -o ${reads.baseName}_nanoplot -f svg --N50 \
  --title "${reads.baseName} sample" --plots hex dot pauvre kde ;

  # Checking Quality
  nanoQC -o ${reads.baseName}_nanoQC ${reads} ;

  # Generate Statistics Summary
  NanoStat --fastq ${reads} -n ${reads.baseName}.txt --outdir ${reads.baseName}_stats ;
  """
}
