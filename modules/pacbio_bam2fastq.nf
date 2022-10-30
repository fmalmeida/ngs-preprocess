process BAM2FASTQ {
  publishDir "${params.output}/preprocessing_outputs/pacbio/bam2fastq", mode: 'copy'
  tag "${meta.id}"
  label 'process_medium'

  input:
  tuple val(meta), path(subreads)
  file barcodes
  
  output:
  tuple val(meta), path("*.fq.gz"), emit: reads
  path "*", emit: all
  
  when:
  !(subreads =~ /input.*/)

  script:
  
  if (params.pacbio_barcodes)
  """
  # index bam
  pbindex ${subreads} ;

  # split bams
  lima ${meta.design} \\
      --num-threads ${task.cpus} \\
      --split-named ${subreads} \\
      ${barcodes} ${meta.id}_demuxed.bam

  # split fastqs
  for input_demux_bam in \$(ls ${meta.id}_demuxed*.bam) ; do
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
  bam2fastq -o ${meta.id} -u ${subreads}

  # fix read extensions and gzip
  for i in *.fastq ; do mv \$i \${i%%.fastq}.fq; gzip \${i%%.fastq}.fq; done
  """
}
