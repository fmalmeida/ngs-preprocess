#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Generic Pipeline for preprocessing ngs data
 */

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/
WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    USEFUL FUNCTION
========================================================================================
*/

// Load / filter data from SRA_FETCH module --> for long reads.
def sra_lreads_filter(in_channel, technology) {
  if (technology == 'nanopore') {
      first_filter_ch = in_channel.filter{ it[1] =~ /nanopore/ }
  } else if (technology == 'pacbio') {
      first_filter_ch = in_channel.filter{ it[1] =~ /pacbio/   }
  }
  result_ch = 
  first_filter_ch
  .map{
      def meta            = [:]
      meta.id             = it[0]
      meta.longreads_type = technology

      [ meta, it[3] ]
  }

  return result_ch
}

/*
========================================================================================
    LOAD INPUT DATA
========================================================================================
*/

//
// SRA
//
sra_ids_ch = (params.sra_ids) ? Channel.value( file(params.sra_ids).readLines() ) : Channel.empty() // if .empty(), workflow does not start

//
// PACBIO DATA
//
subreads_bam      = (params.pacbio_bam)      ? Channel.fromPath(params.pacbio_bam)      : Channel.value('') // if .empty(), workflow does not start
subreads_h5       = (params.pacbio_h5)       ? Channel.fromPath(params.pacbio_h5)       : Channel.value('') // if .empty(), workflow does not start
subreads_barcodes = (params.pacbio_barcodes) ? Channel.fromPath(params.pacbio_barcodes) : Channel.value('') // if .empty(), workflow does not start

//
// NANOPORE DATA
//
sequencing_summary = (params.nanopore_sequencing_summary) ? Channel.fromPath(params.nanopore_sequencing_summary) : Channel.value('') // if .empty(), workflow does not start
if (params.nanopore_fastq) {
  nanopore_fastq =
  Channel.fromPath(params.nanopore_fastq)
  .map{ 
    def meta = [:]
    meta.id = it.getBaseName() - ".fastq.gz" - ".fastq" - ".fq.gz" - ".fq"
    meta.longreads_type = 'nanopore'

    [ meta, it ]
  }
} else {
  nanopore_fastq = Channel.value('') // if .empty(), workflow does not start
}

//
// SHORT READS DATA
//
if (params.shortreads) {
  if (params.shortreads_type == 'paired') {
    shortreads_ch = 
    Channel.fromFilePairs( params.shortreads )
    .map{
      def meta = [:]
      meta.id  = it[0]
      meta.shortreads_type = 'paired'

      [ meta, it[1] ]
    }
  } 
  
  if (params.shortreads_type == 'single') {
    shortreads_ch = 
    Channel.fromPath( params.shortreads )
    .map{
      def meta = [:]
      meta.id  = it.getBaseName() - ".fastq.gz" - ".fastq" - ".fq.gz" - ".fq"
      meta.shortreads_type = 'single'

      [ meta, it ]
    }
  }
} else {
  shortreads_ch = Channel.value('') // if .empty(), workflow does not start
}

/*
========================================================================================
    LOAD WORKFLOWS / MODULES
========================================================================================
*/
include { SRA_FETCH } from './workflows/sra_fetch'
include { NANOPORE  } from './workflows/nanopore'
include { PACBIO    } from './workflows/pacbio'
include { ILLUMINA  } from './workflows/illumina'

/*
========================================================================================
    DEFINE MAIN WORKFLOW
========================================================================================
*/
workflow {

  /*
   * User gives a list of SRA IDs
   */
  SRA_FETCH( sra_ids_ch.flatten() )

  /*
   * User has nanopore longreads
   */
  nanopore_fastq = 
  nanopore_fastq
  .filter{ it != '' }
  .mix( sra_lreads_filter( SRA_FETCH.out.fastqs, 'nanopore' ) )

  NANOPORE( nanopore_fastq, sequencing_summary )

  /*
   * User has pacbio subreads
   */
  pacbio_fastq = sra_lreads_filter( SRA_FETCH.out.fastqs, 'pacbio' )
  PACBIO(
    subreads_h5, 
    subreads_bam, 
    subreads_barcodes,
    pacbio_fastq
  )
  
  /*
   * User has illumina reads
   */
  shortreads_ch = 
  shortreads_ch
  .filter{ it != '' }
  .mix( 
    SRA_FETCH.out.fastqs
    .filter{ it[1] =~ /illumina|bgiseq/ }
    .map{ 
      def meta             = [:]
      meta.id              = it[0]
      meta.shortreads_type = it[2]

      [ meta, it[3] ]
    }
  )
  ILLUMINA( shortreads_ch )

  /*
   * Collect all the generated results as a samplesheet for MpGAP
   */
  
  // get string of final dir
  def final_outdir = file("${params.output}").toUriString()
  final_outdir = "${final_outdir}/final_output"

  // start samplesheet channel and feed it
  ch_mpgap_samplesheet = Channel.value('samplesheet:')
  ch_mpgap_samplesheet.concat(

    // short reads data
    ILLUMINA.out.reads
    .map{ meta, subdir, reads ->
      def reads_list = meta.shortreads_type == 'paired' ? "\s\s\s\s\s\s- ${final_outdir}/${subdir}/${reads[0].getName()}\n\s\s\s\s\s\s- ${final_outdir}/${subdir}/${reads[1].getName()}" : "\s\s\s\s\s\s- ${final_outdir}/${subdir}/${reads.getName()}"
      def final_string = "\n\s\s- id: ${meta.id}"
      final_string = final_string + "\n\s\s\s\sillumina:\n"
      final_string = final_string + reads_list

      final_string
    }

    // nanopore data
  )
  .collectFile( name: 'mpgap_samplesheet.yml', storeDir: params.output, sort: false, cache: false )

}

workflow.onComplete {
    println ""
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
    println "Thank you for using fmalmeida/ngs-preprocess pipeline!"
}


