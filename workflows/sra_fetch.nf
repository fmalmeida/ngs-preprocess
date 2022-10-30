/*
 * Include modules
 */
include { GET_FASTQ    } from '../modules/get_fastq'
include { GET_METADATA } from '../modules/get_metadata'

// def workflow
workflow SRA_FETCH {
  take:
    sra_ids
  
  main:

    // get fastqs
    GET_FASTQ( sra_ids )

    // get sample metadata
    headers = "Run,LibraryName,LibraryStrategy,LibrarySelection,LibraryLayout,SRAStudy,BioProject,Sample,SampleName,Platform"
    GET_METADATA( GET_FASTQ.out.sra_fastq.flatten(), headers )

    // parse results
    GET_METADATA.out.sra_metadata
    .splitCsv( header: true )
    .map{ runInfo ->
      
      if (runInfo.LibraryLayout.toLowerCase() =~ /PAIRED/) {
        regex = "_{1,2}*"
      } else {
        regex = "*"
      }
      
      fastqs_ch = file( "${runInfo.fastq_dir}/${regex}" )

      [ runInfo.Run, runInfo.Platform.toLowerCase(), runInfo.LibraryLayout.toLowerCase(), fastqs_ch ]

    }
    .set{ fastqs_ch }

    emit:
    fastqs = fastqs_ch

}