def exampleMessage() {
   log.info """

   Example Usages:

      Illumina paired end reads. Since it will always be a pattern match, example "illumina/SRR9847694_{1,2}.fastq.gz",
      it MUST ALWAYS be double quoted as the example below.

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir outputs/illumina_paired --shortreads \
"illumina/SRR9847694_{1,2}.fastq.gz" --shortreads_type "paired" --lighter --lighter_genome_size 4600000 \
--clip_r1 5 --three_prime_clip_r1 5 --clip_r2 5 --three_prime_clip_r2 5 --quality_trim 30 --flash


      Illumina single end reads. Multiple files at once, using fixed number of bases to be trimmed
      If multiple unpaired reads are given as input at once, pattern MUST be double quoted: "SRR9696*.fastq.gz"

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir sample_dataset/outputs/illumina_single \
--shortreads "sample_dataset/illumina/SRR9696*.fastq.gz" --shortreads_type "single" --clip_r1 5 --three_prime_clip_r1 5


      ONT reads (filtering reads by length and quality):

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir sample_dataset/outputs/ont \
--nanopore_fastq sample_dataset/ont/kpneumoniae_25X.fastq --lreads_min_length 500 --lreads_min_quality 10


      Pacbio raw (subreads.bam) reads with nextflow general report (filtering reads by length and quality)

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir sample_dataset/outputs/pacbio --pacbio_get_hifi \
--pacbio_bam sample_dataset/pacbio/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreads.bam \
--lreads_min_length 500 --lreads_min_quality 15


      Pacbio raw (legacy .bas.h5 to subreads.bam) reads

./nextflow run fmalmeida/ngs-preprocess --pacbio_h5 E01_1/Analysis_Results/ \
--outdir E01_1/Analysis_Results/preprocessed --threads 3

   """.stripIndent()
}
