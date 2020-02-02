process pacbio_bam2fastq {
  publishDir "${params.outdir}/longreads/pacbio", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Extracting FASTQ from pacbio subreads.bam files"

  input:
    file reads
  output:
    file "*.fastq"

  script:
  id = (reads.getBaseName() - ".bam")
  param = (params.pacbio_is_barcoded) ? "-u --split-barcodes ${reads}" : "-u ${reads}"
  """
  pbindex ${reads} ;
  source activate pbtools ;
  bam2fastq -o ${id} ${param}
  """
}
