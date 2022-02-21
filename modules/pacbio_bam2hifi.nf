process BAM2HIFI {
  publishDir "${params.output}/preprocessing_outputs/pacbio/bam2hifi", mode: 'copy'
  tag "${id}"
  label 'process_medium'

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

  # compute ccs
  ccs \\
      --num-threads ${task.cpus} \\
      ${subreads} \\
      ${id}.ccs.bam

  # split bams
  lima \\
      ${design} \\
      --num-threads ${task.cpus} \\
      --split-named \\
      ${id}.ccs.bam \\
      ${barcodes} \\
      ${id}_demuxed.bam

  # split fastqs
  for input_demux_bam in \$(ls ${id}_demuxed*.bam) ; do
    prefix=\${input_demux_bam%%.bam} ;
    bam2fastq -o \$prefix -u \$input_demux_bam ;
  done

  # fix read extensions and gzip
  for i in *.fastq ; do mv \$i \${i%%.fastq}.fq; gzip \${i%%.fastq}.fq; done
  """
  else
  """
  # index bam
  pbindex ${subreads} ;

  # compute ccs
  ccs --num-threads ${task.cpus} ${subreads} ${id}.ccs.bam ;

  # convert to fastq
  bam2fastq -o ${id}.ccs -u ${id}.ccs.bam
  
  # fix read extensions and gzip
  for i in *.fastq ; do mv \$i \${i%%.fastq}.fq; gzip \${i%%.fastq}.fq; done
  """
}
