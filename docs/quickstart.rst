.. _quickstart:

Quickstart
**********

Overview
--------

We will use few test cases for evaluating the pipeline's commands and workflow.
We've made available two datasets:

* `Dataset 1 <https://drive.google.com/file/d/1xm3R97HXcfhsjSyhoTnvbA8HDK4Ij8ws/view?usp=sharing>`_.

    * Oxford Nanopore data (FAST5 and FASTQ);
    * Illumina paired end reads;

* `Dataset 2 <https://github.com/PacificBiosciences/DevNet/wiki/E-coli-K12-MG1655-Resequencing>`_

    * Pacbio data (subreads.bam);
    * Illumina paired end reads;

Getting the data
================

Users can download the data with the command below.

Add to your ``bashrc`` or ``bash_aliases``
""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

  function gdrive_download () {
  CONFIRM=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate \
  "https://docs.google.com/uc?export=download&id=$1" -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')
  wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$CONFIRM&id=$1" -O $2
  rm -rf /tmp/cookies.txt
  }

.. tip::

  The function is used as: ``gdrive_download [gdrive id] [output name]``

Then, you can download the datasets as follows:

* **Dataset 1**

    * ``gdrive_download 1xm3R97HXcfhsjSyhoTnvbA8HDK4Ij8ws dataset_1.tar.gz``

* **Dataset 2**

    * ``wget http://files.pacb.com/datasets/primary-analysis/e-coli-k12/1.3.0/e-coli-k12-mg1655-raw-reads-1.3.0.tgz``

Preprocessing the data
----------------------

Dataset 1
=========

After downloaded the dataset shall be available as ``dataset_1`` directory. The data can be
preprocessed using the following command:

.. note::

  These parameters can be used via configuration file

Running the pipeline
""""""""""""""""""""

.. code-block:: bash

  # Running for both Illumina and nanopore data
  nextflow run fmalmeida/ngs-preprocess --shortreads "dataset_1/illumina/read_pair_{1,2}.fastq" \
  --shortreads_type "paired" --quality_trim 30 --flash_execute --threads 3 \
  --nanopore_fastq "dataset_1/ont/ont_reads.fastq" --outdir "dataset_1/preprocessed"

Outputs will be at ``dataset_1/preprocessed``

Dataset 2
=========

After downloaded and decompressed the dataset shall be available as ``e-coli-k12-mg1655-raw-reads-1.3.0`` directory. The data can be
preprocessed using the following command:

.. note::

  These parameters can be used via configuration file

Running the pipeline
""""""""""""""""""""

.. code-block:: bash

  # Running for both Illumina and pacbio data
  nextflow run fmalmeida/ngs-preprocess --pacbio_h5Path e-coli-k12-mg1655-raw-reads-1.3.0/2590338/0006/Analysis_Results/*.bas.h5 \
  --outdir e-coli-k12-mg1655-raw-reads-1.3.0/2590338/0006/preprocessed --threads 3

Outputs will be at ``e-coli-k12-mg1655-raw-reads-1.3.0/2590338/0006/preprocessed``

Afterwards
----------

Now you can used these datasets to, for example, assemble a genome.

For this, check `MpGAP <https://mpgap.readthedocs.io/en/latest/index.html>`_ out that we've
developed for assembling reads from Illumina, Pacbio and Oxford Nanopore sequencing platforms.
