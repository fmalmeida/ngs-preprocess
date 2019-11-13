.. _index:

ngs-preprocess A generic pipeline for pre-processing Illumina, Pacbio and ONT data
==================================================================================

`ngs-preprocess <https://github.com/fmalmeida/ngs-preprocess>`_ is a pipeline developed with `Nextflow <https://www.nextflow.io/docs/latest/index.html>`_
and `Docker <https://www.docker.com/>`_. It was designed to provide an easy-to-use framework for performing the pre-process of data from multiple
sequencing platforms (Illumina, Pacbio and Oxford Nanopore). It wraps up the following tools:

* `FastQC <https://github.com/s-andrews/FastQC>`_
* `TrimGalore <https://github.com/FelixKrueger/TrimGalore>`_
* `FLASH <https://ccb.jhu.edu/software/FLASH/>`_
* `Lighter <https://github.com/mourisl/Lighter>`_
* `Porechop <https://github.com/rrwick/Porechop>`_
* `pbh5tools <https://github.com/PacificBiosciences/pbh5tools/blob/master/doc/index.rst>`_
* `bam2fastx <https://github.com/PacificBiosciences/bam2fastx>`_
* `NanoPack <https://github.com/wdecoster/nanopack>`_


.. toctree::
   :hidden:

   installation
   manual
   config
   examples

Support Contact
===============
Whenever a doubt arise feel free to contact me at almeidafmarques@gmail.com
