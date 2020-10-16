.. _index:

ngs-preprocess pipeline
***********************

`ngs-preprocess <https://github.com/fmalmeida/ngs-preprocess>`_ is a pipeline developed with `Nextflow <https://www.nextflow.io/docs/latest/index.html>`_
and `Docker <https://www.docker.com/>`_. It was designed to provide an easy-to-use framework for performing the pre-process of data from multiple
sequencing platforms (Illumina, Pacbio and Oxford Nanopore). It wraps up the following tools:

* `FastQC <https://github.com/s-andrews/FastQC>`_
* `TrimGalore <https://github.com/FelixKrueger/TrimGalore>`_
* `FLASH <https://ccb.jhu.edu/software/FLASH/>`_
* `Lighter <https://github.com/mourisl/Lighter>`_
* `Porechop <https://github.com/rrwick/Porechop>`_
* `bax2bam <https://github.com/PacificBiosciences/bax2bam>`_
* `bam2fastx <https://github.com/PacificBiosciences/bam2fastx>`_
* `lima <https://github.com/PacificBiosciences/barcoding>`_
* `pacbio ccs <https://ccs.how/>`_
* `NanoPack <https://github.com/wdecoster/nanopack>`_
* `pycoQC <https://github.com/a-slide/pycoQC), [bax2bam](https://github.com/PacificBiosciences/bax2bam>`_


.. toctree::
   :hidden:

   installation
   quickstart
   manual
   config
   examples

Support Contact
===============
Feel free to contact me at almeidafmarques@gmail.com
