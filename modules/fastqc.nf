process fastqc {
  publishDir "${params.outdir}/shortreads/before_trimming", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Evaluating short reads with FastQC"

    input:
      file reads
      val threads
    output:
      file "fastqc_${id}/*_fastqc.{zip,html}"
    script:
      if (params.shortreads_type == 'paired') {
        id = (reads[1].getBaseName() - "_1")
        param = "-q ${reads[1]} ${reads[2]}"
      }
      else if (params.shortreads_type == 'single') {
        param = "-q ${reads}"
        id = reads.getBaseName()
      }
    """
    mkdir fastqc_${id} ;
    fastqc -t $threads -o fastqc_${id} $param
    """
}
