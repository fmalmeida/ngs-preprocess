/*
 * Include modules
 */
include { bam2fastq } from '../modules/pacbio_bam2fastq.nf'
include { bam2hifi } from '../modules/pacbio_bam2hifi.nf'
include { h52bam } from '../modules/pacbio_h52bam.nf'
include { nanopack } from '../modules/nanopack.nf'
include { filter } from '../modules/lreads_filter.nf'

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
      h52bam(h5bas)
      h5_bams = h52bam.out
    }

    // has subreads in bam
    // User wants to get hifi?
    reads = Channel.empty()
    if (params.pacbio_get_hifi) {
      bam2hifi(subreads.mix(h5_bams), barcodes)
      reads = bam2hifi.out[0]
    } else {
      bam2fastq(subreads.mix(h5_bams), barcodes)
      reads = bam2fastq.out[0]
    }

    // QC on fastq
    nanopack(reads)

    // filter reads
    if (params.lreads_min_length || params.lreads_min_quality) {
      filter(reads)
    }

}