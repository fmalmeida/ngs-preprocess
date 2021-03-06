.. _quickstart:

Quickstart
**********

Overview
========

As an use case, we will use 30X of one of the *Escherichia coli* sequencing data (Biosample: `SAMN10819847 <https://www.ncbi.nlm.nih.gov/biosample/10819847>`_)
that is available from a recent study that compared the use of different long read technologies in hybrid assembly of 137 bacterial genomes [`4 <https://doi.org/10.1099/mgen.0.000294>`_].

Get the data
============

We have made this subsampled dataset available in `Figshare <https://figshare.com/articles/dataset/Illumina_pacbio_and_ont_sequencing_reads/14036585>`_.

.. code-block:: bash

  # Download data from figshare
  wget -O reads.zip https://ndownloader.figshare.com/articles/14036585/versions/4

  # Unzip
  unzip reads.zip


Now we have the necessary data to perform the quickstart.

.. note::

  The pipeline will always use the fastq file name as prefix for sub-folders and output files. For instance, if users use a fastq file named SRR7128258.fastq the output files and directories will have the string "SRR7128258" in it.

.. tip::

  Remember, the pipeline can always be executed with a config file. In fact, the best way to execute these pipelines is by using a configuration file. With a proper configuration, users can easily run the pipeline.

Preprocessing the data
======================

Outputs will be at ``preprocessed_reads``.

.. warning::

  **Remember:** the pipeline does not concatenate the reads. Whenever you use a pattern such as \* the pipeline will process each pair separately.

.. code-block:: bash

  # Running for both Illumina and nanopore data
  nextflow run fmalmeida/ngs-preprocess \
    --outdir preprocessed_reads \
    --threads 4 \
    --shortreads "SRR8482585_30X_{1,2}.fastq.gz" \
    --shortreads_type paired \
    --lighter_execute \
    --lighter_genomeSize 4m \
    --flash_execute \
    --nanopore_fastq "SRX5299443_30X.fastq.gz" \
    --lreads_min_length 1000 \
    --lreads_min_quality 10

.. note::

  These parameters can be used via configuration file. See :ref:`config`.

Afterwards
==========

Now you can used these datasets to, for example, assemble and annotate a genome. For this, check out the `MpGAP <https://mpgap.readthedocs.io/en/latest/index.html>`_ and `Bacannot <https://bacannot.readthedocs.io/en/latest/index.html>`_ pipelines that we've developed for such tasks.
