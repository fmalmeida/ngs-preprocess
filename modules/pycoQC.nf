process pycoQC {
  publishDir "${params.outdir}/longreads/pycoQC_quality_check", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Checking sequencing statistics with pycoQC"

  input:
    file summary
  output:
    file "pycoQC_report.html"

  script:
    """
    source activate pycoQC ;
    pycoQC -f ${summary} -o pycoQC_report.html --filter_calibration \
    --min_pass_len 200 --min_pass_qual 8 --min_barcode_percent 10
    """
}
