.. _examples:

CLI usage Examples
******************

.. warning::

  Remember: the pipeline does not concatenate the reads. Whenever you use a pattern
  such as \* the pipeline will process each pair separately.

Illumina paired end reads.
""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir outputs/illumina_paired \
      --shortreads "illumina/SRR9847694_{1,2}.fastq.gz" --shortreads_type "paired" \
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
  --nanopore_fastq sample_dataset/ont/kpneumoniae_25X.fastq

Pacbio basecalled (.fastq) reads with nextflow general report
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/pacbio_from_fastq --run_longreads_pipeline --lreads_type pacbio
  --longReads sample_dataset/pacbio/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreads.fastq -with-report

Pacbio raw (subreads.bam) reads
"""""""""""""""""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir sample_dataset/outputs/pacbio \
  --pacbio_bamPath sample_dataset/pacbio/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreads.bam -with-report

Pacbio raw (legacy .bas.h5 to subreads.bam) reads
"""""""""""""""""""""""""""""""""""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --pacbio_h5Path E01_1/Analysis_Results/ \
  --outdir E01_1/Analysis_Results/preprocessed --threads 3

.. note::

  Pacbio bas.h5 file and its related bax.h5 files MUST be in the same directory


Running with a configuration file
"""""""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/ngs-preprocess -c nextflow.config
