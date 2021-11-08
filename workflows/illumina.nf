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

    // run fastqc
    fastqc(preads, sreads)

    // trim reads
    trimgalore(preads, sreads)

    // get trimmed reads
    paired_reads   = (params.shortreads_type == 'paired') ? trimgalore.out[0] : Channel.value(['', '', ''])
    unpaired_reads = (params.shortreads_type == 'single') ? trimgalore.out[1] : Channel.value('')

    // run lighter
    if (params.lighter) {
      lighter(paired_reads, unpaired_reads)
    }

    // run flash
    flash_input = (params.lighter) ? lighter.out[0] : paired_reads
    if (params.flash) {
      flash(flash_input)
    }
}