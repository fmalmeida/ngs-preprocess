process pacbio_bam2fastq {
  publishDir "${params.outdir}/longreads/pacbio", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Extracting FASTQ from pacbio subreads.bam files"

  input:
    file subreads
  output:
    file "*.fastq"

  script:
  id = (subreads.getBaseName() - ".bam")
  param = (params.pacbio_is_barcoded) ? "-u --split-barcodes ${subreads}" : "-u ${subreads}"
  """
  source activate pbtools ;
  pbindex ${subreads} ;
  bam2fastq -o ${id} ${param}
  """
}
