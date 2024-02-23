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
    PORECHOP(reads)
    FILTER(PORECHOP.out.fastqs)
    NANOPACK(PORECHOP.out.fastqs)
  
  emit:
    reads = FILTER.out.reads
}