//
// This file holds several functions specific to the workflow/ngs-preprocess.nf in the fmalmeida/ngs-preprocess pipeline
//

class WorkflowNGSPreprocess {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (!params.get_config && !params.help) {
            if (!params.shortreads && !params.nanopore_fastq && !params.pacbio_bam && !params.pacbio_h5) {
                log.error "Please provide at least on input file! E.g. --shortreads / --nanopore_fastq / --pacbio_bam / --pacbio_h5. Please check the manual at: https://ngs-preprocess.readthedocs.io/en/latest/manual.html#"
                System.exit(1)
            }

            if (params.shortreads) {
                if (params.shortreads_type != 'paired' || params.shortreads_type != 'single') {
                    log.error "Major error: When processing shortreads, params.shortreads_type must be set to either: 'paired' or 'single'. Please, review your parameters."
                    System.exit(1)
                }
            }
        }
    }

}
