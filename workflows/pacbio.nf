/*
 * Include modules
 */
include { BAM2FASTQ } from '../modules/pacbio_bam2fastq.nf'
include { BAM2HIFI  } from '../modules/pacbio_bam2hifi.nf'
include { H52BAM    } from '../modules/pacbio_h52bam.nf'
include { NANOPACK  } from '../modules/nanopack.nf'
include { FILTER    } from '../modules/lreads_filter.nf'

// def workflow
workflow PACBIO {
  take:
    h5bas
    subreads
    barcodes
  main:

    // has H5 data
    h5_bams = Channel.empty()
    if (params.pacbio_h5) {
      H52BAM(h5bas)
      h5_bams = H52BAM.out
    }

    // has subreads in bam
    // User wants to get hifi?
    reads = Channel.empty()
    if (params.pacbio_get_hifi) {
      BAM2HIFI(subreads.mix(h5_bams), barcodes)
      reads = BAM2HIFI.out[0]
    } else {
      BAM2FASTQ(subreads.mix(h5_bams), barcodes)
      reads = BAM2FASTQ.out[0]
    }

    // QC on fastq
    NANOPACK(reads)

    // filter reads
    FILTER(reads)

}