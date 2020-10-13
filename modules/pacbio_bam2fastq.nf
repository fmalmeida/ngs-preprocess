process pacbio_bam2fastq {
  publishDir "${params.outdir}/longreads/pacbio", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Extracting FASTQ from pacbio subreads.bam files"

  input:
    file subreads
    file barcodes
  output:
    file "*.fastq"
    file "*.bam"
    file "*"

  script:
  id = (subreads.getBaseName() - ".bam")
  design = (params.pacbio_barcode_design.toLowerCase() != 'same' && params.pacbio_barcode_design.toLowerCase() != 'different') ? '' : params.pacbio_barcode_design.toLowerCase()
  if (params.pacbio_barcodes)
  """
  source activate pbtools ;
  pbindex ${subreads} ;
  # Split bams
  lima --${design} --num-threads ${params.threads} --split-named ${subreads} ${barcodes} demuxed.bam

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
  bam2fastq -o ${id} -u ${subreads}
  """
}
