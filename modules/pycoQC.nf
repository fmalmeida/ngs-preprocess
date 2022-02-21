process PYCOQC {
  publishDir "${params.output}/preprocessing_outputs/nanopore/QC", mode: 'copy'
  tag "${id}"
  label 'process_low'

  input:
  file summary
  
  output:
  path "${id}_pycoQC.html"

  when:
  !(summary =~ /input.*/)

  script:
  id = summary.getBaseName()
  """
  # run pycoQC
  pycoQC \\
      --summary_file ${summary} \\
      --html_outfile ${id}_pycoQC.html \\
      --filter_calibration \\
      --filter_duplicated \\
      --min_pass_qual 8
  """
}
