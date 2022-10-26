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
    sra_reads

  main:

    // has H5 data
    h5_bams = Channel.empty()
    if (params.pacbio_h5) {
      H52BAM(h5bas)
      h5_bams = H52BAM.out.subreads
    }

    // has subreads in bam
    // User wants to get hifi?
    reads = Channel.empty()
    parsed_subreads = 
    subreads.mix(h5_bams).filter{ it != '' }
    .map{
      def meta    = [:]
      meta.id     = it.getBaseName() - ".bam" - ".subreads"
      meta.design = (params.pacbio_barcode_design.toLowerCase() != 'same' && params.pacbio_barcode_design.toLowerCase() != 'different') ? '' : '--' + params.pacbio_barcode_design.toLowerCase()

      [ meta, it ]
    }
    .view()
    if (params.pacbio_get_hifi) {
      BAM2HIFI(parsed_subreads, barcodes)
      reads = reads.mix( BAM2HIFI.out.reads )
    } else {
      BAM2FASTQ(parsed_subreads, barcodes)
      reads = reads.mix( BAM2FASTQ.out.reads )
    }
    reads = reads.mix(sra_reads)

    // QC on fastq
    NANOPACK(reads)

    // filter reads
    FILTER(reads)

}