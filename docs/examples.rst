.. _examples:

CLI usage Examples
******************

.. warning::

  Remember: the pipeline does not concatenate the reads. Whenever you use a pattern
  such as \* the pipeline will process each pair separately.

.. tip::

  The parameters `--use_tower` and `--tower_token` allows the user to launch the pipeline via `nextflow tower <https://tower.nf>`_ in order to visualize its execution.

Illumina paired end reads.
""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir outputs/illumina_paired \
      --shortreads "illumina/SRR*_{1,2}.fastq.gz" --shortreads_type "paired" \
      --lighter_execute --lighter_genomeSize 4600000 --clip_r1 5 --three_prime_clip_r1 5 \
      --clip_r2 5 --three_prime_clip_r2 5 --quality_trim 30 --flash_execute

.. note::

  Since it will always be a pattern match, example "illumina/SRR9847694_{1,2}.fastq.gz", it MUST ALWAYS be double quoted as the example below.

Illumina single end reads.
""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir sample_dataset/outputs/illumina_single \
      --shortreads "sample_dataset/illumina/SRR9696*.fastq.gz" --shortreads_type "single" --clip_r1 5 --three_prime_clip_r1 5

.. note::

  Multiple files at once, using fixed number of bases to be trimmed. If multiple unpaired reads are given as input at once, pattern MUST be double quoted: "SRR9696*.fastq.gz"

ONT reads (fastq)
"""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir sample_dataset/outputs/ont \
  --nanopore_fastq sample_dataset/ont/kpneumoniae_25X.fastq --lreads_min_length 1000

.. note::

  The parameter ``--lreads_min_length`` applies a min. read length threshold to filter the reads.

Pacbio raw (subreads.bam) reads
"""""""""""""""""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir sample_dataset/outputs/pacbio \
  --pacbio_bamPath sample_dataset/pacbio/m140905_*.subreads.bam -with-report --pacbio_get_hifi

.. note::

  The parameter `--pacbio_get_hifi` will make the pipeline **try** to produce the high fidelity pacbio ccs reads.

Pacbio raw (legacy .bas.h5 to subreads.bam) reads
"""""""""""""""""""""""""""""""""""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --pacbio_h5Path E01_1/Analysis_Results/ \
  --outdir E01_1/Analysis_Results/preprocessed --threads 3

.. note::

  Pacbio bas.h5 file and its related bax.h5 files MUST be in the same directory

Running with a nf-core interactive graphical interface
""""""""""""""""""""""""""""""""""""""""""""""""""""""

::

      ./nf-core launch fmalmeida/ngs-preprocess


Running with a configuration file
"""""""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/ngs-preprocess -c nextflow.config
