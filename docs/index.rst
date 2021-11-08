.. _index:

.. image:: lOGO_3_transparente.png
  :width: 250
  :align: left
  :alt: Laboratory logo

----

NGS-preprocess
**************

`NGS-preprocess <https://github.com/fmalmeida/ngs-preprocess>`_ is a pipeline developed with `Nextflow <https://www.nextflow.io/docs/latest/index.html>`_
and `Docker <https://www.docker.com/>`_. It was designed to provide an easy-to-use framework for preprocessing sequencing reads from Illumina, Pacbio and Oxford Nanopore platforms.

It wraps up the following tools:

.. list-table::
   :widths: 10 60 40
   :header-rows: 1

   * - Software
     - Analysis step
     - Source

   * - TrimGalore
     - Trimming and quality control of Illumina reads
     - https://github.com/FelixKrueger/TrimGalore

   * - Lighter
     - Illumina reads error correction
     - https://github.com/mourisl/Lighter

   * - FLASH
     - Illumina paired end read merger
     - https://ccb.jhu.edu/software/FLASH/

   * - FastQC
     - Illumina reads QC
     - https://github.com/s-andrews/FastQC

   * - Porechop
     - ONT reads trimming and demultiplexing
     - https://github.com/rrwick/Porechop

   * - pycoQC
     - ONT reads QC
     - https://github.com/tleonardi/pycoQC

   * - NanoPack
     - Long reads QC and filter
     - https://github.com/wdecoster/nanopack

   * - bax2bam
     - Convert PacBio bax files to bam
     - https://github.com/PacificBiosciences/bax2bam

   * - bam2fastx
     - Extract reads from PacBio bam files
     - https://github.com/PacificBiosciences/bam2fastx

   * - lima
     - PacBio reads demultiplexing
     - https://github.com/PacificBiosciences/barcoding

   * - pacbio ccs
     - Generate PacBio Highly Accurate Single-Molecule Consensus Reads
     - https://ccs.how/


.. toctree::
   :hidden:

   installation
   profiles
   quickstart
   manual
   config
   examples

Support Contact
***************

Feel free to contact me at almeidafmarques@gmail.com
