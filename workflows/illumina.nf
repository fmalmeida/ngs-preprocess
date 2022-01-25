/*
 * Include modules
 */
include { FASTP } from '../modules/fastp.nf'

// def workflow
workflow ILLUMINA {
  take:
    preads
    sreads
  main:

    // run fastqc
    FASTP(preads, sreads)

}