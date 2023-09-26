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
    
  * - sra-tools & entrez-direct
    - Interaction with SRA database for fetching fastqs and metadata
    - https://anaconda.org/bioconda/entrez-direct ; https://github.com/ncbi/sra-tools

  * - Fastp
    -  tool designed to provide fast all-in-one preprocessing for FastQ files
    - https://github.com/OpenGene/fastp

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
    - https://anaconda.org/bioconda/bax2bam

  * - bam2fastx
    - Extract reads from PacBio bam files
    - https://github.com/PacificBiosciences/pbtk#bam2fastx

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

Citation
********

In order to cite this pipeline, please refer to:

.. code-block:: none
  :class: wrap

      Almeida FMd, Campos TAd and Pappas Jr GJ. Scalable and versatile container-based pipelines for de novo genome assembly and bacterial annotation. [version 1; peer review: awaiting peer review]. F1000Research 2023, 12:1205 (https://doi.org/10.12688/f1000research.139488.1)

Additionally, archived versions of the pipeline are also found in `Zenodo <https://doi.org/10.5281/zenodo.3451405>`_.

Support Contact
***************

Feel free to contact me at almeidafmarques@gmail.com
