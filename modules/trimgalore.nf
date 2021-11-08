process trimgalore {
  publishDir "${params.outdir}/shortreads/${id}/after_trimming", mode: 'copy',
      saveAs: {filename ->
  // This line saves the files with specific sufixes in specific folders
          if (filename.indexOf("_fastqc") > 0) "trim_galore/quality_check/$filename"
          else if (filename.indexOf("trimming_report.txt") > 0) "trim_galore/quality_check/$filename"
          else if (filename.indexOf(".fq.gz") > 0) "trim_galore/$filename"
          else null
      }
  tag "Executing TrimGalore"

  input:
  tuple val(id), file(read1), file(read2)
  file(sreads)

  output:
  tuple val('trimgalore'), file("${id}_1.fq.gz"), file("${id}_2.fq.gz") optional true
  file "*.fq.gz" optional true
  file "*trimming_report.txt"
  file "*_fastqc.{zip,html}"

  when:
  (!(read1 =~ /input.*/) && !(read2 =~ /input.*/)) || !(sreads =~ /input.*/)

  script:
  // Loads Optional Parameters
  c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
  c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
  tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
  tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''

  if (params.shortreads_type == 'paired') {
    param = "--paired $c_r1 $c_r2 $tpc_r1 $tpc_r2 ${read1} ${read2}"
    rename = "mv *_val_1.fq.gz ${id}_1.fq.gz ; mv *_val_2.fq.gz ${id}_2.fq.gz"
  }
  else if (params.shortreads_type == 'single') {
    param = "$c_r1 $tpc_r1 ${sreads}"
    id = reads.getBaseName()
    rename = ''
  }
  """
  # run trim_galore
  trim_galore -q ${params.quality_trim} --fastqc --gzip $param -j ${params.threads} ;
  ${rename}
  """
}
