/*
 * Include modules
 */
include { FASTP } from '../modules/fastp.nf'

// def workflow
workflow ILLUMINA {
  take:
    shortreads
  
  main:

    // run fastqc
    FASTP( shortreads )

}