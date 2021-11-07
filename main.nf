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
//include { paramsCheck } from './nf_functions/paramsCheck.nf'
include { configMessage; illuminaMessage; ontMessage; pacbioMessage; logMessage } from './nf_functions/logMessages.nf'

/*
 * Check parameters
 */
//paramsCheck()
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
params.get_full_config = false
params.get_illumina_config = false
params.get_ont_config = false
params.get_pacbio_config = false

if (params.get_full_config) {
  new File("ngs-preprocess.config").write(new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/nextflow.config").getText())
  configMessage()
  exit 0
}

if (params.get_illumina_config) {
  new File("illumina_data.config").write(new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/configuration_example/illumina_data.config").getText())
  illuminaMessage()
  exit 0
}

if (params.get_ont_config) {
  new File("ont_data.config").write(new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/configuration_example/ont_data.config").getText())
  ontMessage()
  exit 0
}

if (params.get_pacbio_config) {
  new File("pacbio_data.config").write(new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/configuration_example/pacbio_data.config").getText())
  pacbioMessage()
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
params.lighter_execute = false
params.lighter_kmer = 21
params.lighter_genomeSize = false
params.lighter_alpha = false
params.flash_execute = false

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



include { fastqc } from './modules/fastqc.nf' params(outdir: params.outdir,
  shortreads_type: params.shortreads_type)

include { trimgalore } from './modules/trimgalore.nf' params(outdir: params.outdir,
  shortreads_type: params.shortreads_type, clip_r1: params.clip_r1,
  clip_r2: params.clip_r2, three_prime_clip_r1: params.three_prime_clip_r1,
  three_prime_clip_r2: params.three_prime_clip_r2, quality_trim: params.quality_trim)

include { lighter } from './modules/lighter.nf' params(outdir: params.outdir,
  lighter_kmer: params.lighter_kmer, lighter_alpha: params.lighter_alpha,
  shortreads_type: params.shortreads_type, lighter_genomeSize: params.lighter_genomeSize)

include { flash } from './modules/flash.nf' params(outdir: params.outdir)

/*
 * Define custom workflows
 */

// workflow shortreads_nf {
//   take:
//     preads
//     sreads
//     threads
//   main:
//     fastqc(preads, sreads, threads)
//     trimgalore(preads, sreads, threads)

//     // Paired
//     if (params.shortreads_type == 'paired') {
//       if (params.lighter_execute && params.flash_execute) {
//         lighter(trimgalore.out[0], threads)
//         flash(lighter.out[0], threads)
//       } else if (params.lighter_execute && !params.flash_execute) {
//         lighter(trimgalore.out[0], threads)
//       } else if (!params.lighter_execute && params.flash_execute) {
//         flash(trimgalore.out[0], threads)
//       }
//     }

//     // Single
//     if (params.shortreads_type == 'single') {
//       if (params.lighter_execute) {
//         lighter(trimgalore.out[1].flatten(), threads)
//       }
//     }
// }

/*
 * Define main workflow
 */
workflow {
  /*
   * User has nanopore longreads
   */
  sequencing_summary = (params.nanopore_sequencing_summary) ? Channel.fromPath(params.nanopore_sequencing_summary) : ''
  nanopore_fastq = (params.nanopore_fastq) ? Channel.fromPath(params.nanopore_fastq) : ''
  if (params.nanopore_fastq) {
    NANOPORE(nanopore_fastq, sequencing_summary)
  }

  /*
   * User has pacbio subreads
   */
  subreads_bam = (params.pacbio_bam) ? Channel.fromPath(params.pacbio_bam) : ''
  subreads_h5  = (params.pacbio_h5)  ? Channel.fromPath(params.pacbio_h5)  : ''
  subreads_barcodes = (params.pacbio_barcodes) ? Channel.fromPath(params.pacbio_barcodes) : ''
  if (params.pacbio_h5 || params.pacbio_bam) {
    PACBIO(subreads_h5, subreads_bam, subreads_barcodes)
  }
  
  // /*
  //  * User has short paired end reads
  //  */
  // if (params.shortreads && params.shortreads_type == 'paired') {
  //   shortreads_nf(Channel.fromFilePairs(params.shortreads, flat: true, size: 2),
  //                 Channel.value(''), params.threads)
  // }

  // /*
  //  * User has short single end reads
  //  */
  // if (params.shortreads && params.shortreads_type == 'single') {
  //   shortreads_nf(Channel.value(['', '', '']), Channel.fromPath(params.shortreads), params.threads)
  // }
}

workflow.onComplete {
    println ""
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
    println "Thank you for using fmalmeida/ngs-preprocess pipeline!"
}
