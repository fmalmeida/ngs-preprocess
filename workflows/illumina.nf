/*
 * Include modules
 */
include { fastqc } from '../modules/fastqc.nf'
include { trimgalore } from '../modules/trimgalore.nf'
include { lighter } from '../modules/lighter.nf'
include { flash } from '../modules/flash.nf'

// def workflow
workflow ILLUMINA {
  take:
    preads
    sreads
  main:
    fastqc(preads, sreads)
    trimgalore(preads, sreads)

    // Paired
    if (params.shortreads_type == 'paired') {
      if (params.lighter && params.flash) {
        lighter(trimgalore.out[0])
        flash(lighter.out[0])
      } else if (params.lighter && !params.flash) {
        lighter(trimgalore.out[0])
      } else if (!params.lighter && params.flash) {
        flash(trimgalore.out[0])
      }
    }

    // Single
    if (params.shortreads_type == 'single') {
      if (params.lighter) {
        lighter(trimgalore.out[1].flatten())
      }
    }
}