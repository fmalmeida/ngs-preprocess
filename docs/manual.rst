.. _manual:

Manual
******

Input
=====

    * path to fastq files containing sequencing reads (Illumina or Nanopore)
    * path to Pacbio .bam or .h5 files containing raw data

.. note::

   Users must **never** use hard or symbolic links. This will make nextflow fail.
   When setting the parameters, please **always** give full path to a hard file,
   not to a link. This will prevent file access fail.

.. warning::

  Remember: the pipeline does not concatenate the reads. Whenever you use a pattern
  such as \* the pipeline will process each pair separately.

.. tip::

  The parameters `--use_tower` and `--tower_token` allows the user to launch the pipeline via
  `nextflow tower <https://tower.nf>`_ in order to visualize its execution.


Usage example
=============

::

   nextflow run fmalmeida/ngs-preprocess [OPTIONS]

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--use_tower``
     - N
     - False
     - Triggers the pipeline to be launched via nextflow tower

   * - ``--tower_token``
     - Y (if ``--use_tower``)
     - NA
     - Your nextflow tower token. Used to launch the pipeline in your nextflow tower account.

   * - ``--outdir``
     - Y
     - output
     - Name of directory to store output values

   * - ``--threads``
     - N
     - 2
     - Number of threads to use

   * - ``--shortreads``
     - Y
     - NA
     - String Pattern to find short reads. Example: "SRR6307304_{1,2}.fastq"

   * - ``--shortreads_type``
     - Y
     - NA
     - (single | paired). Tells wheter input is unpaired or paired end.

   * - ``--clip_r1``
     - N
     - 0
     - Number of bases to always remove from 5' of read pair 1 or from unpaired read.

   * - ``--clip_r2``
     - N
     - 0
     - Number of bases to always remove from 5' of read pair 2.

   * - ``--three_prime_clip_r1``
     - N
     - 0
     - Number of bases to always remove from 3' of read pair 1 or from unpaired read

   * - ``--three_prime_clip_r2``
     - N
     - 0
     - Number of bases to always remove from 3' of read pair 2.

   * - ``--quality_trim``
     - N
     - 20
     - Phred quality threshold for trimming.

   * - ``--lighter_execute``
     - N
     - False
     - Tells wheter to run or not Lighter correction tool

   * - ``--lighter_kmer``
     - N
     - 21
     - Lighter k-mer to use in correction step.

   * - ``--lighter_genomeSize``
     - Y (If ``--lighter_execute``)
     - NA
     - Approximate genome size

   * - ``--lighter_alpha``
     - N
     - NA
     - Lighter sample rate alpha parameter. If empty, Lighter will automatically calculate its value.

   * - ``--flash_execute``
     - N
     - False
     - If set, FLASH will be executed to merge paired end reads

   * - ``--lreads_min_length``
     - N
     - NA
     - Length min. threshold for filtering long reads (ONT or Pacbio).

   * - ``--lreads_min_quality``
     - N
     - NA
     - Quality min. threshold for filtering long reads (ONT or Pacbio).

   * - ``--nanopore_fastq``
     - Y
     - NA
     - Sets path to nanopore fastq files. Pre-processes basecalled long reads.

   * - ``--nanopore_is_barcoded``
     - N
     - False
     - Tells wheter your data (Nanopore or Pacbio) is barcoded or not. It will split barcodes into single files. Users with legacy pacbio data need to first produce a new barcoded_subreads.bam file.

   * - ``--nanopore_sequencing_summary``
     - N
     - NA
     - Path to nanopore 'sequencing_summary.txt'. Using this will make the pipeline render a sequencing statistics report using pycoQC

   * - ``--pacbio_bamPath``
     - N
     - NA
     - Path to Pacbio subreads.bam. Only used if user wants to basecall subreads.bam to FASTQ. Always keep subreads.bam and its relative subreads.bam.pbi files in the same directory

   * - ``--pacbio_h5Path``
     - N
     - NA
     - Path to directory containing legacy bas.h5 data file (1 per directory). It will be used to extract reads in FASTQ file. All its related files (e.g. bax.h5 files) must be in the same directory

   * - ``--pacbio_barcodes``
     - N
     - False
     - Path to xml/fasta file containing barcode information. It will split barcodes into single files.

   * - ``--pacbio_barcode_design``
     - N
     - same
     - Select the combination of barcodes for demultiplexing. Options: same, different, any.

   * - ``--pacbio_get_hifi``
     - N
     - False
     - Whether or not to try to compute CCS reads


All this parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Examples
--------

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
