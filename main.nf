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
include { paramsCheck } from './nf_functions/paramsCheck.nf'
include { configMessage; illuminaMessage; ontMessage; pacbioMessage; logMessage } from './nf_functions/logMessages.nf'

/*
 * Check parameters
 */
paramsCheck()
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
params.outdir = 'output'
params.threads = 2

/*
 * Parameters for short reads
 */
params.shortreads = ''
params.shortreads_type = 'paired' //paired or single
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0
params.quality_trim = 20
params.lighter = false
params.lighter_kmer = 21
params.lighter_genome_size = false
params.lighter_alpha = false
params.flash = false

/*
 * Parameters for longreads filtering
 */
params.lreads_min_quality = false
params.lreads_min_length  = false

/*
 * Parameters for nanopore longreads
 */
params.nanopore_fastq = ''
params.nanopore_is_barcoded = false
params.nanopore_sequencing_summary = ''

/*
 * Parameters for pacbio longreads
 */
params.pacbio_bam = ''
params.pacbio_h5 = ''
params.pacbio_barcodes = ''
params.pacbio_barcode_design = 'same'
params.pacbio_get_hifi = false

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
