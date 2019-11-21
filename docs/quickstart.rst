.. _quickstart:

Quickstart
==========

Example dataset
---------------

We provide users a few test cases for evaluating the pipeline's commands and workflow.
We've made available two datasets:

* `Dataset 1 <https://drive.google.com/file/d/1OJImuNgNQo_Wxbi3QnErcQXzhwuEPQNL/view?usp=sharing>`_.

    * Oxford Nanopore data (FAST5 and FASTQ);
    * Illumina paired end reads;

* `Dataset 2 <>`_

    * Pacbio data (FASTA and subreads.bam);
    * Illumina paired end reads;

Getting the data
----------------

Users can directly download data through the link given or by this command line method below.

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

    * ``gdrive_download 1OJImuNgNQo_Wxbi3QnErcQXzhwuEPQNL kp_ont_dataset.tar.gz``

* **Dataset 2**

Preprocessing the data
----------------------

Illumina reads
~~~~~~~~~~~~~~
