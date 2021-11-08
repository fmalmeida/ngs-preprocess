process lighter {
  publishDir "${params.outdir}/shortreads/${id}/after_trimming", mode: 'copy',
  saveAs: {filename ->
  // This line saves the files with specific sufixes in specific folders
           if (filename.indexOf(".cor{.fq.gz, .fq}") > 0) "lighter_corrected/$filename"
           else "lighter_corrected/$filename"}
  tag "Executing Ligther (read correction) step"

  input:
  tuple val(pair_id), file(read1), file(read2)
  tuple val(unpaired_id), (sreads)

  output:
  tuple val(id), file("*_1.cor.fq.gz"), file("*_2.cor.fq.gz") optional true
  tuple val(id), file("*.cor.fq.gz") optional true
  file 'fastqc_after_correction'

  when:
  (!(read1 =~ /input.*/) && !(read2 =~ /input.*/)) || !(sreads =~ /input.*/)

  script:
  // Check if alpha is given
  alpha_param = (params.lighter_alpha) ? "-k ${params.lighter_kmer} ${params.lighter_genome_size} ${params.lighter_alpha}" : "-K ${params.lighter_kmer} ${params.lighter_genome_size}"
  // Check reads library
  if (params.shortreads_type == 'paired') {
      param = "-r ${read1} -r ${read2}"
      quality = "*_1.cor.fq.gz *_2.cor.fq.gz"
      id = pair_id
  } else if (params.shortreads_type == 'single') {
      param = "-r ${sreads}"
      quality = "*.cor{.fq.gz, .fq}"
      id = unpaired_id
  }
  """
  # run lighter
  lighter -t ${params.threads} ${param} ${alpha_param};
  mkdir fastqc_after_correction ;
  fastqc -t ${params.threads} -o fastqc_after_correction -q ${quality}
  """
}
