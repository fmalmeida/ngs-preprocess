def configMessage() {
  log.info """

  ngs-preprocess.config file saved in working directory
  After configuration, run:
  nextflow run fmalmeida/ngs-preprocess -c ./ngs-preprocess.config
  Nice code

  """.stripIndent()
}

/*
 * Define log message
 */
def logMessage() {
  log.info "================================="
  log.info "fmalmeida/ngs-preprocess pipeline"
  log.info "================================="
  def summary = [:]
  summary['Output dir']   = params.output
  if(workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Current home']   = "$HOME"
  //summary['Current user']   = "$USER"
  summary['Current path']   = "$PWD"
  log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
  log.info "==================================="
}
