/*
 * Include modules
 */
include { PORECHOP     } from '../modules/porechop.nf'
include { PORECHOP_ABI } from '../modules/porechop_abi.nf'
include { NANOPACK     } from '../modules/nanopack.nf'
include { FILTER       } from '../modules/lreads_filter.nf'
include { PYCOQC       } from '../modules/pycoQC.nf'

// def workflow
workflow NANOPORE {
  take:
    reads
    fast5
  
  main:
    if (params.nanopore_sequencing_summary) {
      PYCOQC(fast5)
    }
    if (params.use_porechop_abi) {
      PORECHOP_ABI(reads)
      ch_nano = PORECHOP_ABI.out.fastqs
    } else {
      PORECHOP(reads)
      ch_nano = PORECHOP.out.fastqs
    }
    FILTER(ch_nano)
    NANOPACK(ch_nano)
  
  emit:
    reads = FILTER.out.reads
}