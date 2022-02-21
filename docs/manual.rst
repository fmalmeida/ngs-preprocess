.. _manual:

Manual
******

Input
=====

* path to fastq files containing sequencing reads (Illumina or Nanopore)
* path to Pacbio .bam or .h5 files containing raw data

.. warning::

  Users must **never** use hard or symbolic links. This will make nextflow fail.

  Whenever using REGEX for a pattern match, for example "illumina/SRR9847694_{1,2}.fastq.gz" or "illumina/SRR*.fastq.gz", it MUST ALWAYS be inside double quotes.

  **Remember:** the pipeline does not concatenate the reads. Whenever you use a pattern such as \* the pipeline will process each read (or pair) that match this pattern separately.


Usage example
=============

.. code-block:: bash

   nextflow run fmalmeida/ngs-preprocess [--help] [OPTIONS]

Output directory
================

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--output``
     - Y
     - output
     - Name of directory to store output values


Max job request
===============

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--max_cpus``
     - N
     - 4
     - Max number of threads to use in parallel
   
   * - ``--max_memory``
     - N
     - 6.GB
     - Max amount of memory to be used by pipeline
   
   * - ``--max_time``
     - N
     - 40.h
     - Max time for a job to run


Short reads (Illumina)
======================

.. list-table::
   :widths: 30 10 10 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--shortreads``
     - Y
     - NA
     - String Pattern to find short reads. Example: "SRR6307304_{1,2}.fastq"

   * - ``--shortreads_type``
     - Y
     - NA
     - (single | paired). Tells whether input is unpaired or paired end.

   * - ``--fastp_average_quality``
     - N
     - 20
     - Fastp will filter out reads with mean quality less than this.

   * - ``--fastp_correct_pairs``
     - N
     - False
     - If set, tells Fastp to try to correct paired end reads. Only works for paired end reads.

   * - ``--fastp_merge_pairs``
     - N
     - False
     - If set, tells Fastp to try to merge read pairs.

   * - ``--fastp_additional_parameters``
     - N
     - False
     - Pass on any additional parameter to Fastp. The tool's parameters are described in their `manual <https://github.com/OpenGene/fastp>`_


Long reads (Pacbio or Nanopore)
===============================

.. list-table::
   :widths: 30 10 10 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--lreads_min_length``
     - N
     - 500
     - Length min. threshold for filtering long reads (ONT or Pacbio).

   * - ``--lreads_min_quality``
     - N
     - 5
     - Quality min. threshold for filtering long reads (ONT or Pacbio).

   * - ``--nanopore_fastq``
     - Y
     - NA
     - Sets path to nanopore fastq files. Pre-processes basecalled long reads.

   * - ``--nanopore_is_barcoded``
     - N
     - False
     - Tells whether your data (Nanopore or Pacbio) is barcoded or not. It will split barcodes into single files. Users with legacy pacbio data need to first produce a new barcoded_subreads.bam file.

   * - ``--nanopore_sequencing_summary``
     - N
     - NA
     - Path to nanopore 'sequencing_summary.txt'. Using this will make the pipeline render a sequencing statistics report using pycoQC. pycoQC reports will be saved using the files basename, so please, use meaningful basename, such as: sample1.txt, sample2.txt, etc. Preferentially, using the same basename as the fastq.

   * - ``--pacbio_bam``
     - N
     - NA
     - Path to Pacbio subreads.bam. Only used if user wants to basecall subreads.bam to FASTQ. Always keep subreads.bam and its relative subreads.bam.pbi files in the same directory

   * - ``--pacbio_h5``
     - N
     - NA
     - Path to directory containing legacy bas.h5 data file (1 per directory). It will be used to extract reads in FASTQ file. All its related files (e.g. bax.h5 files) must be in the same directory

   * - ``--pacbio_barcodes``
     - N
     - False
     - Path to xml/fasta file containing barcode information. It will split barcodes into single files. Will be used for all pacbio inputs, h5 or bam.

   * - ``--pacbio_barcode_design``
     - N
     - same
     - Select the combination of barcodes for demultiplexing. Options: same, different, any.

   * - ``--pacbio_get_hifi``
     - N
     - False
     - Whether or not to try to compute CCS reads. Will be used for all pacbio inputs, h5 or bam.


All this parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Examples
========

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
