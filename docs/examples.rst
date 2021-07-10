.. _examples:

CLI usage Examples
******************

.. warning::

  Remember: the pipeline does not concatenate the reads. Whenever you use a pattern
  such as \* the pipeline will process each pair separately.

Illumina paired end reads.
""""""""""""""""""""""""""

.. code-block:: bash

      ./nextflow run fmalmeida/ngs-preprocess \
        --threads 3 \
        --outdir illumina_paired \
        --shortreads "path-to/SRR*_{1,2}.fastq.gz" \
        --shortreads_type "paired" \
        --flash_execute

.. note::

  Since ``--shortreads`` will always be a pattern match, example "illumina/SRR9847694_{1,2}.fastq.gz", it MUST ALWAYS be double quoted as the example below.

.. note::

  ``--flash_execute`` triggers paired end reads merge with FLASH.

Illumina single end reads.
""""""""""""""""""""""""""

.. code-block:: bash

      ./nextflow run fmalmeida/ngs-preprocess \
        --threads 3 \
        --outdir illumina_single \
        --shortreads "path-to/SRR9696*.fastq.gz" \
        --shortreads_type "single" \
        --clip_r1 5 --three_prime_clip_r1 5

.. note::

  Multiple files at once, using fixed number of bases to be trimmed (``--clip_r1`` and ``--three_prime_clip_r1``). If multiple unpaired reads are given as input at once, pattern MUST be double quoted: "SRR9696*.fastq.gz"

ONT reads (fastq)
"""""""""""""""""

.. code-block:: bash

  ./nextflow run fmalmeida/ngs-preprocess \
    --threads 3 \
    --outdir ONT \
    --nanopore_fastq "path-to/ont_reads.fastq" \
    --lreads_min_length 1000

.. note::

  The parameter ``--lreads_min_length`` applies a min. read length threshold to filter the reads.

Pacbio raw (subreads.bam) reads
"""""""""""""""""""""""""""""""

.. code-block:: bash

  ./nextflow run fmalmeida/ngs-preprocess \
    --threads 3 \
    --outdir pacbio_subreads \
    --pacbio_bamPath "path-to/m140905_*.subreads.bam" \
    --pacbio_get_hifi \
    -with-report

.. note::

  The parameter ``--pacbio_get_hifi`` will make the pipeline **try** to produce the high fidelity pacbio ccs reads.

.. note::

  ``-with-report`` will generate nextflow execution reports.

Pacbio raw (legacy .bas.h5 to subreads.bam) reads
"""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

  ./nextflow run fmalmeida/ngs-preprocess \
    --pacbio_h5Path E01_1/Analysis_Results/ \
    --outdir E01_1/Analysis_Results/preprocessed \
    --threads 3

.. note::

  This example refers to the SMRT Cell data files available at: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly. The path ``E01_1/Analysis_Results/`` is the directory where the legacy \*.bas.h5 and \*.bax.h5 files are located. The pipeline will load the bas files available in the directory.

.. note::

  Pacbio bas.h5 file and its related bax.h5 files MUST be in the same directory

Running with a nf-core interactive graphical interface
""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: bash

      ./nf-core launch fmalmeida/ngs-preprocess


Running with a configuration file
"""""""""""""""""""""""""""""""""

.. code-block:: bash

      ./nextflow run fmalmeida/ngs-preprocess -c nextflow.config
