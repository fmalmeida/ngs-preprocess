process trimgalore {
    publishDir "${params.outdir}/illumina/after_trimming", mode: 'copy',
        saveAs: {filename ->
    // This line saves the files with specific sufixes in specific folders
            if (filename.indexOf("_fastqc") > 0) "quality_check/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "quality_check/$filename"
            else if (filename.indexOf(".fq.gz") > 0) "reads/trimmed/$filename"
            else null
        }
    container 'fmalmeida/ngs-preprocess'
    tag "Executing TrimGalore with paired end reads."

    input:
      file reads
      file threads
    output:
      tuple val('trimgalore'), file("${id}_1.fq.gz"), file("${id}_2.fq.gz") optional true
      file "*.fq.gz" optional true
      file "*trimming_report.txt"
      file "*_fastqc.{zip,html}"

    script:
      // Loads Optional Parameters
      c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
      c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
      tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
      tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''

      if (params.shortreads_type == 'paired') {
        id = (reads[1].getBaseName() - "_1")
        param = "--paired $c_r1 $c_r2 $tpc_r1 $tpc_r2 ${reads[1]} ${reads[2]}"
        rename = "mv *_val_1.fq.gz ${id}_1.fq.gz ; mv *_val_2.fq.gz ${id}_2.fq.gz"
      }
      else if (params.shortreads_type == 'single') {
        param = "$c_r1 $tpc_r1 ${reads}"
        id = reads.getBaseName()
        rename = ''
      }

    """
    trim_galore -q ${params.quality_trim} --fastqc --gzip $param ;
    ${rename}
    """
}
