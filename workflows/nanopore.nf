/*
 * Include modules
 */
include { porechop } from '../modules/porechop.nf'
include { nanopack } from '../modules/nanopack.nf'
include { filter } from '../modules/lreads_filter.nf'
include { pycoQC } from '../modules/pycoQC.nf'

// def workflow
workflow NANOPORE {
  take:
    reads
    fast5
  main:
    if (params.nanopore_sequencing_summary) {
      pycoQC(fast5)
    }
    if (params.nanopore_fastq){
      porechop(reads)
      // barcoded reads
      barcoded = (params.nanopore_is_barcoded) ? porechop.out[1] : Channel.value('')
      not_barcoded = !(params.nanopore_is_barcoded) ? porechop.out[0] : Channel.value('')
      nanopack(not_barcoded.mix(barcoded))
      if (params.lreads_min_length || params.lreads_min_quality) {
        filter(not_barcoded.mix(barcoded))
      }
    }
}