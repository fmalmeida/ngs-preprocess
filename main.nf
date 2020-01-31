#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
 * Define help messages
 */
def helpMessage() {
   log.info """

   Usage:
   nextflow run fmalmeida/ngs-preprocess [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]

   Comments:
   This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would
   cause the command to be huge.

   Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make
   parameterization easier and more readable.

   Creating a configuration file:
   nextflow run fmalmeida/ngs-preprocess [--get_illumina_config] [--get_ont_config] [--get_pacbio_config]

   Show command line examples:
   nextflow run fmalmeida/ngs-preprocess --examples

   Execution Reports:
   nextflow run fmalmeida/ngs-preprocess [ -c nextflow.config ] -with-report
   nextflow run fmalmeida/ngs-preprocess [ -c nextflow.config ] -with-trace
   nextflow run fmalmeida/ngs-preprocess [ -c nextflow.config ] -with-timeline

   OBS: These reports can also be enabled through the configuration file.

   OPTIONS:

            General Parameters - Mandatory

    --outDir <string>                      Output directory name
    --threads <int>                        Number of threads to use
    --run_shortreads_pipeline              Selects preprocess pipeline of Illumina short reads
    --run_longreads_pipeline               Selects preprocess pipeline of ONT or Pacbio long reads

            Short Reads Parameters - Mandatory if --run_shortreads_pipeline is used

    --shortreads <string>                  String Pattern to find short reads. Example: SRR6307304_{1,2}.fastq
    --reads_size <int>                     Tells wheter input is unpaired or paired end. 1 is unpaired. 2 is paired
    --clip_r1 <int>                        Number of bases to always remove from 5' of read pair 1 or from unpaired read. [Default: 0]
    --clip_r2 <int>                        Number of bases to always remove from 5' of read pair 2. [Default: 0]
    --three_prime_clip_r1 <int>            Number of bases to always remove from 3' of read pair 1 or from unpaired read. [Default: 0]
    --three_prime_clip_r2 <int>            Number of bases to always remove from 3' of read pair 2. [Default: 0]
    --quality_trim <int>                   Phred quality threshold for trimming. [Default: 20]
    --lighter_execute                      Tells wheter to run or not Lighter correction tool
    --lighter_kmer <int>                   Lighter k-mer to use in correction step. [Default: 21]
    --lighter_genomeSize <int>             Approximate genome size
    --lighter_alpha <float>                Lighter sample rate alpha parameter. Rule of thumb: (7/C) where C is coverage.
                                           If not set, Lighter will automatically calculate the best value
    --flash_execute                        If set, PEAR will be executed to merge paired end reads

            Long Reads Parameters - Mandatory if --run_shortreads_pipeline is used

    --lreads_type <string>                 Tells wheter input is from pacbio or ONT. Possibilities: pacbio | nanopore
    --longReads <string>                   Path to pacbio or ONT basecalled reads. If reads are already basecalled, the pipeline
                                           will only check sequencing quality and statistics with NanoPack
    --lreads_is_barcoded                   Tells wheter your data is barcoded or not. It will split barcodes into single files.
                                           Users with legacy pacbio data need to first produce a new barcoded_subreads.bam file.

            For PacificBiosciences Data - Mandatory if --lreads_type pacbio

    --pacbio_bamPath <string>              Path to Pacbio subreads.bam. Only used if user wants to basecall subreads.bam to FASTQ.
                                           Always keep subreads.bam and its relative subreads.bam.pbi files in the same directory
    --pacbio_h5Path <string>               Path to legacy *.bas.h5 data. It will be used to extract reads in FASTQ file.
                                           All related *bas.h5 and *bax.h5 files MUST be in the SAME dir.

   """.stripIndent()
}

def exampleMessage() {
   log.info """

   Example Usages:

      Illumina paired end reads. Since it will always be a pattern match, example "illumina/SRR9847694_{1,2}.fastq.gz",
      it MUST ALWAYS be double quoted as the example below.

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir outputs/illumina_paired --run_shortreads_pipeline --shortreads \
"illumina/SRR9847694_{1,2}.fastq.gz" --reads_size 2 --lighter_execute --lighter_genomeSize 4600000 --clip_r1 5 --three_prime_clip_r1 5 \
--clip_r2 5 --three_prime_clip_r2 5 --quality_trim 30 --flash_execute


      Illumina single end reads. Multiple files at once, using fixed number of bases to be trimmed
      If multiple unpaired reads are given as input at once, pattern MUST be double quoted: "SRR9696*.fastq.gz"

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/illumina_single --run_shortreads_pipeline \
--shortreads "sample_dataset/illumina/SRR9696*.fastq.gz" --reads_size 1 --clip_r1 5 --three_prime_clip_r1 5


      ONT reads:

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/ont --run_longreads_pipeline \
--lreads_type nanopore --longReads sample_dataset/ont/kpneumoniae_25X.fastq --nanopore_prefix kpneumoniae_25X


      Pacbio basecalled (.fastq) reads with nextflow general report

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/pacbio_from_fastq \
--run_longreads_pipeline --lreads_type pacbio \
--longReads sample_dataset/pacbio/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreads.fastq -with-report


      Pacbio raw (subreads.bam) reads

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/pacbio --run_longreads_pipeline \
--lreads_type pacbio --pacbio_bamPath sample_dataset/pacbio/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreads.bam


      Pacbio raw (legacy .bas.h5 to subreads.bam) reads

./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/pacbio --run_longreads_pipeline \
--lreads_type pacbio --pacbio_h5Path sample_dataset/pacbio/m140912_020930_00114_c100702482550000001823141103261590_s1_p0.1.bas.h5

   """.stripIndent()
}

def illuminaMessage() {
  log.info """

  illumina_data.config file saved in working directory
  After configuration, run:
  nextflow run fmalmeida/ngs-preprocess -c ./illumina_data.config
  Nice code

  """.stripIndent()
}

def ontMessage() {
  log.info """

  ont_data.config file saved in working directory
  After configuration, run:
  nextflow run fmalmeida/ngs-preprocess -c ./ont_data.config
  Nice code

  """.stripIndent()
}

def pacbioMessage() {
  log.info """

  pacbio_data.config file saved in working directory
  After configuration, run:
  nextflow run fmalmeida/ngs-preprocess -c ./pacbio_data.config
  Nice code

  """.stripIndent()
}

/*
 * Check if user want some help
 */
params.help = false
if (params.help){
  helpMessage()
  exit 0
}
params.examples = false
if (params.examples){
  exampleMessage()
  exit 0
}

/*
 * Download configuration files if requested
 */
params.get_illumina_config = false
params.get_ont_config = false
params.get_pacbio_config = false
if (params.get_illumina_config) {
  new File("illumina_data.config") << new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/configuration_example/illumina_data.config").getText()
  illuminaMessage()
  exit 0
}

if (params.get_ont_config) {
  new File("ont_data.config") << new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/configuration_example/ont_data.config").getText()
  ontMessage()
  exit 0
}

if (params.get_pacbio_config) {
  new File("pacbio_data.config") << new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/configuration_example/pacbio_data.config").getText()
  pacbioMessage()
  exit 0
}

/*
 * Load general parameters and establish defaults
 */
params.outdir = 'output'
params.threads = 2
params.run_shortreads_pipeline = false
params.run_longreads_pipeline  = false

/*
 * Parameters for short reads
 */
params.shortreads = ''
params.shortreads_type = 'paired' //paired or single
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0
params.quality_trim = 20
params.lighter_execute = false
params.lighter_kmer = 21
params.lighter_genomeSize = false
params.lighter_alpha = false
params.flash_execute = false

/*
 * Parameters for long reads
 */
params.lreads_type = ''
params.pacbio_bamPath = ''
params.pacbio_h5Path = ''
params.lreads_is_barcoded = false
params.longreads = ''
longreads = (params.longreads) ? Channel.fromPath(params.longreads) : ''

/*
 * Define log message
 */
log.info "==================================="
log.info " fmalmeida/ngs-preprocess pipeline "
log.info "==================================="
def summary = [:]
summary['Output dir']   = params.outdir
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 * Include modules
 */
include porechop from './modules/porechop.nf' params(outdir: params.outdir)
include fastqc from './modules/fastqc.nf' params(outdir: params.outdir,
  shortreads_type: params.shortreads_type)
include trimgalore from './modules/trimgalore.nf' params(outdir: params.outdir,
  shortreads_type: params.shortreads_type, clip_r1: params.clip_r1,
  clip_r2: params.clip_r2, three_prime_clip_r1: params.three_prime_clip_r1,
  three_prime_clip_r2: params.three_prime_clip_r2, quality_trim: params.quality_trim)

/*
 * Define custom workflows
 */
workflow porechop_nf {
  get:
    reads
    threads
    barcode
  main:
    porechop(reads, threads, barcode)
}

workflow fastqc_nf {
  get:
    reads
    threads
  main:
    fastqc(reads, threads)
}

workflow trimgalore_nf {
  get:
    reads
    threads
  main:
    trimgalore(reads, threads)
}

/*
 * Define main workflow
 */
workflow {
  /*
   * User has long reads
   */
  if (params.longreads) {
    porechop_nf(longreads, params.threads, params.lreads_is_barcoded)
  }

  /*
   * User has short paired end reads
   */
  if (params.shortreads && params.shortreads_type == 'paired') {
    fastqc_nf(Channel.fromFilePairs(params.shortreads, flat: true, size: 2), params.threads)
    trimgalore_nf(Channel.fromFilePairs(params.shortreads, flat: true, size: 2), params.threads)
  }

  /*
   * User has short single end reads
   */
  if (params.shortreads && params.shortreads_type == 'single') {
    fastqc_nf(Channel.fromPath(params.shortreads), params.threads)
    trimgalore_nf(Channel.fromPath(params.shortreads), params.threads)
  }
}
