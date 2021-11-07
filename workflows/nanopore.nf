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
    porechop(reads)
    nanopack(porechop.out[0].flatten())
    if (params.lreads_min_length || params.lreads_min_quality) {
      filter(porechop.out[0].flatten())
    }
}