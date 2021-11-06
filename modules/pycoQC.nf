process pycoQC {
  publishDir "${params.outdir}/longreads/pycoQC_quality_check", mode: 'copy'
  tag "Checking sequencing statistics with pycoQC"

  input:
    file summary
  output:
    file "pycoQC_report.html"

  script:
    """
    # run pycoQC
    pycoQC --summary_file ${summary} --html_outfile pycoQC_report.html \
    --filter_calibration --filter_duplicated --min_pass_qual 8
    """
}
