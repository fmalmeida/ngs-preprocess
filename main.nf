#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Generic Pipeline for preprocessing ngs data
 */

/*
 * Include functions
 */
include { helpMessage } from './nf_functions/help.nf'
include { exampleMessage } from './nf_functions/examples.nf'
include { configMessage; logMessage } from './nf_functions/logMessages.nf'

/*
 * Check parameters
 */
params.help = false
if (params.help){
  helpMessage()
  exit 0
}
params.examples = false
if (params.examples){
  exampleMessage()
  exit 0
}

/*
 * Download configuration files if requested
 */
params.get_config = false
if (params.get_config) {
  new File("ngs-preprocess.config").write(new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/nextflow.config").getText())
  configMessage()
  exit 0
}

/*
 * Load general parameters and establish defaults
 */
params.output  = 'output'
params.threads = 2

/*
 * Parameters for short reads
 */
params.shortreads                  = ''
params.shortreads_type             = 'paired' // paired or single
params.fastp_average_quality       = 20
params.fastp_merge_pairs           = false
params.fastp_correct_pairs         = false
params.fastp_additional_parameters = ''

/*
 * Parameters for longreads filtering
 */
params.lreads_min_quality = 5
params.lreads_min_length  = 500

/*
 * Parameters for nanopore longreads
 */
params.nanopore_fastq              = ''
params.nanopore_is_barcoded        = false
params.nanopore_sequencing_summary = ''

/*
 * Parameters for pacbio longreads
 */
params.pacbio_bam            = ''
params.pacbio_h5             = ''
params.pacbio_barcodes       = ''
params.pacbio_barcode_design = 'same'
params.pacbio_get_hifi       = false

/*
 * Define log message
 */
logMessage()

/*
 * Load workflows
 */
include { NANOPORE } from './workflows/nanopore.nf'
include { PACBIO   } from './workflows/pacbio.nf'
include { ILLUMINA } from './workflows/illumina.nf'

/*
 * Define main workflow
 */
workflow {
  /*
   * User has nanopore longreads
   */
  sequencing_summary = (params.nanopore_sequencing_summary) ? Channel.fromPath(params.nanopore_sequencing_summary) : Channel.value('')
  nanopore_fastq = (params.nanopore_fastq) ? Channel.fromPath(params.nanopore_fastq) : Channel.value('')
  if (params.nanopore_fastq || params.nanopore_sequencing_summary) {
    NANOPORE(nanopore_fastq, sequencing_summary)
  }

  /*
   * User has pacbio subreads
   */
  subreads_bam = (params.pacbio_bam) ? Channel.fromPath(params.pacbio_bam) : Channel.value('')
  subreads_h5  = (params.pacbio_h5)  ? Channel.fromPath(params.pacbio_h5)  : Channel.value('')
  subreads_barcodes = (params.pacbio_barcodes) ? Channel.fromPath(params.pacbio_barcodes) : Channel.value('')
  if (params.pacbio_h5 || params.pacbio_bam) {
    PACBIO(subreads_h5, subreads_bam, subreads_barcodes)
  }
  
  /*
   * User has illumina reads
   */
  paired_reads = (params.shortreads && params.shortreads_type == 'paired') ? Channel.fromFilePairs(params.shortreads, flat: true, size: 2) : Channel.value(['', '', ''])
  unpaired_reads = (params.shortreads && params.shortreads_type == 'single') ? Channel.fromPath(params.shortreads) : Channel.value('')
  if (params.shortreads && (params.shortreads_type == 'paired' || params.shortreads_type == 'single')) {
    ILLUMINA(paired_reads, unpaired_reads)
  }
}

workflow.onComplete {
    println ""
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
    println "Thank you for using fmalmeida/ngs-preprocess pipeline!"
}
