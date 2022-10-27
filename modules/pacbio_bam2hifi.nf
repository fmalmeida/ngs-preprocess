process BAM2HIFI {
  publishDir "${params.output}/preprocessing_outputs/pacbio/bam2hifi", mode: 'copy'
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

  # compute ccs
  ccs \\
      --num-threads ${task.cpus} \\
      ${subreads} \\
      ${meta.id}.ccs.bam

  # split bams
  lima \\
      ${meta.design} \\
      --num-threads ${task.cpus} \\
      --split-named \\
      ${meta.id}.ccs.bam \\
      ${barcodes} \\
      ${meta.id}_demuxed.bam

  # split fastqs
  for input_demux_bam in \$(ls ${meta.id}_demuxed*.bam) ; do
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
  ccs --num-threads ${task.cpus} ${subreads} ${meta.id}.ccs.bam ;

  # convert to fastq
  bam2fastq -o ${meta.id}.ccs -u ${meta.id}.ccs.bam
  
  # fix read extensions and gzip
  for i in *.fastq ; do mv \$i \${i%%.fastq}.fq; gzip \${i%%.fastq}.fq; done
  """
}
