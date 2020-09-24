process pacbio_bam2hifi {
  publishDir "${params.outdir}/longreads/pacbio/ccs_hifi", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Computing HIFI reads from pacbio subreads.bam files"

  input:
    file subreads
  output:
    file "*.fastq"
    file "${id}.ccs*"

  script:
  id = (subreads.getBaseName() - ".bam")
  param = (params.pacbio_is_barcoded) ? "-u --split-barcodes ${id}.ccs.bam" : "-u ${id}.ccs.bam"
  """
  source activate pbtools ;
  pbindex ${subreads} ;
  ccs --num-threads ${params.threads} ${subreads} ${id}.ccs.bam ;
  bam2fastq -o ${id}.ccs ${param}
  """
}
