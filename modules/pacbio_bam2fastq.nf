process bam2fastq {
  publishDir "${params.outdir}/longreads/${id}", mode: 'copy'
  tag "Extracting FASTQ from pacbio subreads.bam files"

  input:
  file subreads
  file barcodes
  
  output:
  file "*.fastq"
  file "*"
  
  when:
  !(subreads =~ /input.*/)

  script:
  id = (subreads.getBaseName() - ".bam")
  design = (params.pacbio_barcode_design.toLowerCase() != 'same' && params.pacbio_barcode_design.toLowerCase() != 'different') ? '' : '--' + params.pacbio_barcode_design.toLowerCase()
  if (params.pacbio_barcodes)
  """
  # index bam
  pbindex ${subreads} ;
  # Split bams
  lima ${design} --num-threads ${params.threads} --split-named ${subreads} ${barcodes} demuxed.bam

  # Split fastqs
  for input_demux_bam in \$(ls demuxed*.bam) ; do
    prefix=\${input_demux_bam%%.bam} ;
    # convert bam
    bam2fastq -o \$prefix -u \$input_demux_bam ;
  done
  """
  else
  """
  # index bam
  pbindex ${subreads} ;
  # convert bam
  bam2fastq -o ${id} -u ${subreads}
  """
}
