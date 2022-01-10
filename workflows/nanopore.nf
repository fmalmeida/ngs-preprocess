/*
 * Include modules
 */
include { PORECHOP } from '../modules/porechop.nf'
include { NANOPACK } from '../modules/nanopack.nf'
include { FILTER   } from '../modules/lreads_filter.nf'
include { PYCOQC   } from '../modules/pycoQC.nf'

// def workflow
workflow NANOPORE {
  take:
    reads
    fast5
  main:
    if (params.nanopore_sequencing_summary) {
      PYCOQC(fast5)
    }
    if (params.nanopore_fastq){
      PORECHOP(reads)
      // barcoded reads
      barcoded = (params.nanopore_is_barcoded) ? PORECHOP.out[1].transpose() : Channel.value(['', '', ''])
      not_barcoded = !(params.nanopore_is_barcoded) ? PORECHOP.out[0] : Channel.value(['', '', ''])
      NANOPACK(not_barcoded.mix(barcoded))
      FILTER(not_barcoded.mix(barcoded))
    }
}