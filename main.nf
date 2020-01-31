#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
 * Include modules
 */
include porechop from './modules/porechop.nf'

/*
 * Load parameters and establishing defaults
 */
params.lreads_is_barcoded = false
params.threads = 2
longreads = (params.longreads) ? Channel.fromPath(params.longreads) : ''

/*
 * Define custom workflows
 */
workflow porechop_nf {
  get:
    reads
    threads
    barcode
  main:
    porechop(reads, threads, barcode)
}

/*
 * Define main workflow
 */
workflow {
  main:
    porechop_nf(longreads, params.threads, params.lreads_is_barcoded)
}
