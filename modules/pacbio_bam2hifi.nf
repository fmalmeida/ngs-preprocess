process pacbio_bam2hifi {
  publishDir "${params.outdir}/longreads/pacbio/ccs_hifi", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Computing HIFI reads from pacbio subreads.bam files"

  input:
    file subreads
    file barcodes
  output:
  file "*.fastq"
  file "*"

  script:
  id = (subreads.getBaseName() - ".bam")
  design = (params.pacbio_barcode_design.toLowerCase() != 'same' && params.pacbio_barcode_design.toLowerCase() != 'different') ? '' : '--' + params.pacbio_barcode_design.toLowerCase()
  if (params.pacbio_barcodes)
  """
  source activate pbtools ;
  pbindex ${subreads} ;
  ccs --num-threads ${params.threads} ${subreads} ${id}.ccs.bam

  # Split bams
  lima ${design} --num-threads ${params.threads} --split-named ${id}.ccs.bam ${barcodes} demuxed.bam

  # Split fastqs
  for input_demux_bam in \$(ls demuxed*.bam) ; do
    prefix=\${input_demux_bam%%.bam} ;
    bam2fastq -o \$prefix -u \$input_demux_bam ;
  done
  """
  else
  """
  source activate pbtools ;
  pbindex ${subreads} ;
  ccs --num-threads ${params.threads} ${subreads} ${id}.ccs.bam ;
  bam2fastq -o ${id}.ccs -u ${id}.ccs.bam
  """
}
