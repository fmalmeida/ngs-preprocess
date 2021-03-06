.. _config:

Configuration File
******************

To download a configuration file template users just need to run: ``nextflow run fmalmeida/ngs-preprocess [--get_illumina_config] [--get_ont_config] [--get_pacbio_config]``

Using a config file your code is a lot more cleaner and concise: ``nextflow run fmalmeida/ngs-preprocess -c [path-to-config]``

Default configuration:
""""""""""""""""""""""

.. code-block:: groovy

  /*


                                        fmalmeida/ngs-preprocess pipeline configuration file

                      Maintained by Felipe Marques de Almeida
                      Contact: almeidafmarques@outlook.com

                      This file conatains all the parameters that are required by the pipeline.
                      Some of them can be left in blank and other ones need a proper configuration.
                      Every file containing a parameter will contain comments to guide its configuration.




   */

  /*


                      Customizable parameters

                      All parameters are set under the params {} function.

                      Please do not change anything before the = sign.

                      Fields must be configured and set with values after the = sign.

                                  E.g. param1 = 'user_value_1'


   */

  params {
    /*


                      General (mandatory) Parameters


     */

    outdir  = 'output'                                     // Sets output directory
    threads = 2                                            // Sets number of threads to be used
    parallel_jobs = 1                                      // Number of jobs to run in parallel. Be aware that each job (in parallel) can consume
                                                           // N threads (set above). Be sure to carefully check your resources before augmenting
                                                           // this parameter. For example: parallel_jobs = 2 + threads = 5 can consume until 10
                                                           // threads at once.

    /*


                      Parameters for short reads preprocessing


     */

    shortreads = ''                                       // Path to shortreads. Examples: 'SRR6307304_{1,2}.fastq' | 'SRR7128258*'
    shortreads_type = ''                                  // Type of shortreads. Values: paired | single
    clip_r1 = 0                                           // Number of bases to ALWAYS clip from 5' (read 1) end, despite base qualities
    clip_r2 = 0                                           // Number of bases to ALWAYS clip from 5' (read 2) end, despite base qualities
    three_prime_clip_r1 = 0                               // Number of bases to ALWAYS clip from 3' (read 1) end, despite base qualities
    three_prime_clip_r2 = 0                               // Number of bases to ALWAYS clip from 3' (read 2) end, despite base qualities
    quality_trim = 20                                     // Quality threshold for trimming
    lighter_execute = false                               // Tells whether or not to execute lighter correction step
    lighter_kmer = 21                                     // Which k-mer to use in lighter correction. Check Ligther's manual (https://github.com/mourisl/Lighter)
    lighter_genomeSize = 0                                // Tells lighter the expected genome size for correction of reads
    lighter_alpha = ''                                    // Lighter alpha parameter. Rule of thumb: (7/C) where C is coverage.
                                                          // If left blank, Lighter will automatically calculate the best value.
    flash_execute = false                                 // Tells wheter or not to merge paired reads with FLASH

    /*

                        Optional parameter for filtering longreads from pacbio or nanopore.
                        This command uses nanofilt for both technologies. If you prefer to
                        filter reads from any technology with a different program you can
                        left it blank and it you not be executed.

    */
    lreads_min_quality =                                  // If blank, lreads will not be filtered.
    lreads_min_length  =                                  // If blank, lreads will not be filtered.

    /*


                        Parameters for nanopore ONT longreads preprocessing


     */

    nanopore_fastq = ''                                   // Path to nanopore ONT basecalled reads in fastq
    nanopore_is_barcoded = false                          // Tells wheter or not nanopore reads are barcoded
                                                          // It will split barcodes into single files
    nanopore_sequencing_summary = ''                      // Path to nanopore 'sequencing_summary.txt'. Using this will make the pipeline render a
                                                          // sequencing statistics report using pycoQC

    /*


                        Parameters for PacBio longreads preprocessing

                        Use bamPath or h5Path, not both.


     */

    pacbio_bamPath  = ''                                   // Path to PacBio subreads in bam format
    pacbio_h5Path   = ''                                   // Path to directory containing legacy *.bas.h5 data (1 per directory)
    pacbio_barcodes = ''                                   // Path to xml/fasta file containing barcode information. It will split barcodes into single files.
    pacbio_barcode_design = ''                             // By default, only reads with "same" barcodes are given. You can also select reads with only
                                                           // "different" barcodes or any of them. Options: same, different, any
    pacbio_get_hifi = false                                // Whether or not to try to compute CCS reads

  }
