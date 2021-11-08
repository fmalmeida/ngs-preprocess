process fastqc {
  publishDir "${params.outdir}/shortreads/${id}/before_trimming", mode: 'copy'
  tag "Evaluating short reads with FastQC"

    input:
      tuple val(id), file(read1), file(read2)
      file(sreads)

    output:
      file "fastqc_${id}/*_fastqc.{zip,html}"

    script:

      if (params.shortreads_type == 'paired') {
        param = "-q ${read1} ${read2}"
      }
      else if (params.shortreads_type == 'single') {
        param = "-q ${sreads}"
        id = sreads.getBaseName()
      }
      
    """
    # create directory for results
    mkdir fastqc_${id} ;

    # run fastqc
    fastqc -t ${params.threads} -o fastqc_${id} $param
    """
}
