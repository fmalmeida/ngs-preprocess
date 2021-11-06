def configMessage() {
  log.info """

  ngs-preprocess.config file saved in working directory
  After configuration, run:
  nextflow run fmalmeida/ngs-preprocess -c ./ngs-preprocess.config
  Nice code

  """.stripIndent()
}

def illuminaMessage() {
  log.info """

  illumina_data.config file saved in working directory
  After configuration, run:
  nextflow run fmalmeida/ngs-preprocess -c ./illumina_data.config
  Nice code

  """.stripIndent()
}

def ontMessage() {
  log.info """

  ont_data.config file saved in working directory
  After configuration, run:
  nextflow run fmalmeida/ngs-preprocess -c ./ont_data.config
  Nice code

  """.stripIndent()
}

def pacbioMessage() {
  log.info """

  pacbio_data.config file saved in working directory
  After configuration, run:
  nextflow run fmalmeida/ngs-preprocess -c ./pacbio_data.config
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
  summary['Output dir']   = params.outdir
  if(workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Current home']   = "$HOME"
  //summary['Current user']   = "$USER"
  summary['Current path']   = "$PWD"
  log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
  log.info "==================================="
}
