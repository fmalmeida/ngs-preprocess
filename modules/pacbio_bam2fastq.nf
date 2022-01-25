process BAM2FASTQ {
  publishDir "${params.output}/preprocessing_outputs/pacbio/bam2fastq", mode: 'copy'
  tag "${id}"

  input:
  file subreads
  file barcodes
  
  output:
  tuple val(id), file("*.fq.gz"), val('pacbio')
  file "*"
  
  when:
  !(subreads =~ /input.*/)

  script:
  // variables
  id = (subreads.getBaseName() - ".bam" - ".subreads")
  design = (params.pacbio_barcode_design.toLowerCase() != 'same' && params.pacbio_barcode_design.toLowerCase() != 'different') ? '' : '--' + params.pacbio_barcode_design.toLowerCase()

  // script
  if (params.pacbio_barcodes)
  """
  # index bam
  pbindex ${subreads} ;

  # split bams
  lima ${design} \\
      --num-threads ${params.threads} \\
      --split-named ${subreads} \\
      ${barcodes} ${id}_demuxed.bam

  # split fastqs
  for input_demux_bam in \$(ls ${id}_demuxed*.bam) ; do
    prefix=\${input_demux_bam%%.bam} ;
    # convert bam
    bam2fastq -o \$prefix -u \$input_demux_bam ;
  done

  # fix read extensions and gzip
  for i in *.fastq ; do mv \$i \${i%%.fastq}.fq; gzip \${i%%.fastq}.fq; done
  """
  
  else
  """
  # index bam
  pbindex ${subreads} ;

  # convert bam
  bam2fastq -o ${id} -u ${subreads}

  # fix read extensions and gzip
  for i in *.fastq ; do mv \$i \${i%%.fastq}.fq; gzip \${i%%.fastq}.fq; done
  """
}
