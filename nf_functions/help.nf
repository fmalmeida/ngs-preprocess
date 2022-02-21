/*
 * Define help messages
 */
def helpMessage() {
   log.info """

   Usage:
   nextflow run fmalmeida/ngs-preprocess [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]

   Comments:
   This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would cause the command to be huge.

   Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make parameterization easier and more readable.

   Creating a configuration file:
   nextflow run fmalmeida/ngs-preprocess [--get_config]

   Execution Reports:
   nextflow run fmalmeida/ngs-preprocess [OPTIONS] [-with-report] [-with-trace] [-with-timeline]

   OBS: These reports can also be enabled through the configuration file.
   OBS 2: Make sure parameter values are double quoted

   OPTIONS:

            # General Parameters

    --output <string>                              Output directory name. Default: 'output'
    
    --max_cpus <int>                               Max number of threads to use in parallel. Default: 4
    
    --max_memory <string>                          Max amount of memory to use. E.g. '6.GB'. Default: 6.GB
    
    --max_time <string>                            Max time for a process to run. E.g. '40.h'. Default: 40.h


            # Parameters for short reads preprocessing

    --shortreads <string>                          String Pattern to find short reads. Example: SRR6307304_{1,2}.fastq

    --shortreads_type <string>                     Possibilities: single | paired. Tells wheter input is single or paired end.

    --fastp_average_quality <int>                  Fastp will filter out reads with mean quality less than this. [Default: 20]

    --fastp_merge_pairs                            Tells Fastp to try to merge read pairs.

    --fastp_correct_pairs                          Tells Fastp to try to correct paired end reads.

    --fastp_additional_parameters <string>         Pass on any additional parameter to Fastp.


            # Parameters for long reads filtering
            # Works with both nanopore and pacbio

    --lreads_min_length <int>                      If set, the pipeline will filter the longreads by this minimun length. [Default: 500]

    --lreads_min_quality <int>                     If set, the pipeline will filter the longreads by this minimun quality. [Default: 5]


            # Parameters for preprocessing NANOPORE long reads

    --nanopore_fastq <string>                      Path to ONT basecalled reads.

    --nanopore_is_barcoded                         Inform the pipeline that the data is barcoded. It will split barcodes into single files.

    --nanopore_sequencing_summary                  Path to nanopore 'sequencing_summary.txt'. Using this will make the pipeline render a
                                                   sequencing statistics report using pycoQC


            # Parameters for preprocessing PACBIO long reads
            # PACBIO bam files or legacy h5

    --pacbio_bam <string>                          Path to Pacbio subreads.bam. Only used if user wants to basecall subreads.bam to FASTQ.

    --pacbio_h5 <string>                           Path to directory containing legacy *.bas.h5 data (1 per directory). It will be used to
                                                   extract reads in FASTQ file. All its related files (e.g. bax.h5 files) must be in the same directory.

    --pacbio_barcodes                              Path to xml/fasta file containing barcode information. It will split barcodes into single files.

    --pacbio_barcode_design                        By default, only reads with "same" barcodes are given. You can also select reads with only
                                                   "different" barcodes or any of them. Options: same, different, any

    --pacbio_get_hifi                              Also try to use pbccs to compute subreads consensus and produce HIFI reads. ccs combines multiple subreads
                                                   of the same SMRTbell molecule. Therefore, the bam files used as input must already be merged since this tool
                                                   takes one bam (from one movie) at a time. Can be used for the legacy *.bas.h5 since this pipeline
                                                   automatically creates one subreads.bam for each single movies (each *.bas.h5). If the chemistry is incompatible
                                                   with ccs an error will be thrown and you can re-run the pipeline removing this parameter, using '-resume'.

   """.stripIndent()
}
