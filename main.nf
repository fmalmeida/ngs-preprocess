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
params.pacbio_bamPath = ''
params.pacbio_h5Path = ''
params.pacbio_barcodes = ''
params.pacbio_barcode_design = 'same'
params.pacbio_get_hifi = false

/*
 * Define log message
 */
logMessage()

/*
 * Include modules
 */
include { porechop } from './modules/porechop.nf' params(outdir: params.outdir)

include { nanopack; nanopack as nanopack_hifi } from './modules/nanopack.nf' params(outdir: params.outdir)

include { lreads_filter; lreads_filter as lreads_filter_hifi } from './modules/lreads_filter.nf' params(outdir: params.outdir,
  lreads_min_length: params.lreads_min_length, lreads_min_quality: params.lreads_min_quality)

include { pycoQC } from './modules/pycoQC.nf' params(outdir: params.outdir)

include { pacbio_bam2fastq } from './modules/pacbio_bam2fastq.nf' params(outdir: params.outdir,
  pacbio_barcodes: params.pacbio_barcodes, pacbio_barcode_design: params.pacbio_barcode_design,
  threads: params.threads)

include { pacbio_bam2hifi } from './modules/pacbio_bam2hifi.nf' params(outdir: params.outdir,
  pacbio_barcodes: params.pacbio_barcodes, pacbio_barcode_design: params.pacbio_barcode_design,
  threads: params.threads)

include { pacbio_h52bam } from './modules/pacbio_h52bam.nf' params(outdir: params.outdir)

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
workflow nanopore_nf {
  take:
    reads
    threads
    barcode
  main:
    porechop(reads, threads, barcode)
    nanopack(porechop.out[0].flatten(), threads)
    if (params.lreads_min_length || params.lreads_min_quality) {
      lreads_filter(porechop.out[0].flatten())
    }
}

workflow pycoQC_nf {
  take:
    input
  main:
    pycoQC(input)
}

workflow pacbio_bam_nf {
  take:
    subreads
    barcodes
    threads
  main:
    pacbio_bam2fastq(subreads, barcodes)
    nanopack(pacbio_bam2fastq.out[0].flatten(), threads)
    if (params.lreads_min_length || params.lreads_min_quality) {
      lreads_filter(pacbio_bam2fastq.out[0].flatten())
    }

    // User wants to get hifi?
    if (params.pacbio_get_hifi) {
      pacbio_bam2hifi(subreads, barcodes)
      nanopack_hifi(pacbio_bam2hifi.out[0].flatten(), threads)
      if (params.lreads_min_length || params.lreads_min_quality) {
        lreads_filter_hifi(pacbio_bam2hifi.out[0].flatten())
      }
    }
}

workflow pacbio_bas_nf {
  take:
    h5bas
    h5bas_dir
    barcodes
    threads
  main:
    pacbio_h52bam(h5bas, h5bas_dir)
    bams = pacbio_h52bam.out[0]
    pacbio_bam2fastq(bams, barcodes)
    nanopack(pacbio_bam2fastq.out[0].flatten(), threads)
    if (params.lreads_min_length || params.lreads_min_quality) {
      lreads_filter(pacbio_bam2fastq.out[0].flatten())
    }

    // User wants to get hifi?
    if (params.pacbio_get_hifi) {
      pacbio_bam2hifi(bams, barcodes)
      nanopack_hifi(pacbio_bam2hifi.out[0].flatten(), threads)
      if (params.lreads_min_length || params.lreads_min_quality) {
        lreads_filter_hifi(pacbio_bam2hifi.out[0].flatten())
      }
    }
}

workflow shortreads_nf {
  take:
    preads
    sreads
    threads
  main:
    fastqc(preads, sreads, threads)
    trimgalore(preads, sreads, threads)

    // Paired
    if (params.shortreads_type == 'paired') {
      if (params.lighter_execute && params.flash_execute) {
        lighter(trimgalore.out[0], threads)
        flash(lighter.out[0], threads)
      } else if (params.lighter_execute && !params.flash_execute) {
        lighter(trimgalore.out[0], threads)
      } else if (!params.lighter_execute && params.flash_execute) {
        flash(trimgalore.out[0], threads)
      }
    }

    // Single
    if (params.shortreads_type == 'single') {
      if (params.lighter_execute) {
        lighter(trimgalore.out[1].flatten(), threads)
      }
    }
}

/*
 * Define main workflow
 */
workflow {
  /*
   * User has nanopore longreads
   */
  if (params.nanopore_fastq) {
    nanopore_nf(Channel.fromPath(params.nanopore_fastq), params.threads, params.nanopore_is_barcoded)
  }

  /*
   * User wants to render a report with pycoQC
   */
  if (params.nanopore_sequencing_summary) {
    pycoQC_nf(Channel.fromPath(params.nanopore_sequencing_summary))
  }

  /*
   * User has pacbio subreads in bam format
   */
  if (params.pacbio_bamPath) {
    pacbio_bam_nf(Channel.fromPath(params.pacbio_bamPath),
                  (params.pacbio_barcodes) ? Channel.fromPath(params.pacbio_barcodes) : Channel.value(''),
                  params.threads)
  }

  /*
   * User has pacbio subreads in legacy h5 (bas and bax) files
   */
  if (params.pacbio_h5Path) {
    pacbio_bas_nf(Channel.fromPath(params.pacbio_h5Path),
                  Channel.fromPath(params.pacbio_h5Path, type: 'dir'),
                  (params.pacbio_barcodes) ? Channel.fromPath(params.pacbio_barcodes) : Channel.value(''),
                  params.threads)
  }

  /*
   * User has short paired end reads
   */
  if (params.shortreads && params.shortreads_type == 'paired') {
    shortreads_nf(Channel.fromFilePairs(params.shortreads, flat: true, size: 2),
                  Channel.value(''), params.threads)
  }

  /*
   * User has short single end reads
   */
  if (params.shortreads && params.shortreads_type == 'single') {
    shortreads_nf(Channel.value(['', '', '']), Channel.fromPath(params.shortreads), params.threads)
  }
}

workflow.onComplete {
    println ""
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
    println "Thank you for using fmalmeida/ngs-preprocess pipeline!"
}
