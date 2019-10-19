#!/usr/bin/env nextflow

/*
 * A docker-based pipeline to run the preprocess steps for Illumina, Pacbio and Oxford
 * Nanopore Technologies raw data.
 * It is meant to cut addaptors using TrimGalore for Illumina data and porechop
 * for long reads data.
 *
 * For more details read the docs at:
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
    --flash_execute                         If set, PEAR will be executed to merge paired end reads

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

/*

          Display Help Message

*/

params.help = false
 // Show help emssage
 if (params.help){
   helpMessage()
   file('work').deleteDir()
   file('.nextflow').deleteDir()
   exit 0
}

// CLI examples
params.examples = false
 // Show help emssage
 if (params.examples){
   exampleMessage()
   exit 0
}

/*

          Defaults. They are meant to be overwriten.
          Through CLI or configuration file

*/

params.outDir = ''
params.threads = 2
params.run_shortreads_pipeline = false
params.run_longreads_pipeline  = false
params.shortreads = ''
params.reads_size = 2
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
params.lreads_type = ''
params.longReads = ''
params.lreads_is_barcoded = false
params.pacbio_bamPath = ''
params.pacbio_h5Path = ''
params.get_illumina_config = false
params.get_ont_config = false
params.get_pacbio_config = false

/*

          Download configuration file, if necessary.

*/

if (params.get_illumina_config) {
  new File("illumina_data.config") << new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/configuration_example/illumina_data.config").getText()
  println ""
  println "illumina_data.config file saved in working directory"
  println "After configuration, run:"
  println "nextflow run fmalmeida/ngs-preprocess -c ./illumina_data.config"
  println "Nice code!\n"

  exit 0
}

if (params.get_ont_config) {
  new File("ont_data.config") << new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/configuration_example/ont_data.config").getText()
  println ""
  println "ont_data.config file saved in working directory"
  println "After configuration, run:"
  println "nextflow run fmalmeida/ngs-preprocess -c ./ont_data.config"
  println "Nice code!\n"

  exit 0
}

if (params.get_pacbio_config) {
  new File("pacbio_data.config") << new URL ("https://github.com/fmalmeida/ngs-preprocess/raw/master/configuration_example/pacbio_data.config").getText()
  println ""
  println "pacbio_data.config file saved in working directory"
  println "After configuration, run:"
  println "nextflow run fmalmeida/ngs-preprocess -c ./pacbio_data.config"
  println "Nice code!\n"

  exit 0
}

/*

          STEP 0 - Load configuration parameters

*/

// General Parameters
  // Output Directory
outdir = params.outDir

// Short Reads Parameters
  // Loading Paired or Unpaired Short Reads.
if (params.reads_size == 2 && params.run_shortreads_pipeline) {

    // Loading Paired End Short Reads

    Channel.fromFilePairs(params.shortreads, flat: true, size: 2).into { fastqc_paired; trimgalore_paired }

    // Create Empty Channels if input is not Paired End Short Reads.

    Channel.empty().into { fastqc_single; trimgalore_single }
  } else if (params.reads_size == 1 && params.run_shortreads_pipeline) {

    // Loading Unpaired Short Reads

    Channel.fromPath(params.shortreads).into { fastqc_single; trimgalore_single }

    // Create Empty Channels if input is not Unpaired Short Reads.

    Channel.empty().into {fastqc_paired; trimgalore_paired }
  } else {
    Channel.empty().into { fastqc_paired; trimgalore_paired; fastqc_single; trimgalore_single }
  }

// Long Reads Parameters

/*
      Oxford Nanopore Reads - ONLY

      Loading ONT fastq files
*/

lreads = (params.longReads && params.run_longreads_pipeline && params.lreads_type == 'nanopore') ?
         Channel.fromPath(params.longReads) : Channel.empty()

/*

      PacificBiosciences data - ONLY

      Loading pacbio subreads.bam files

*/

if (params.pacbio_bamPath && params.run_longreads_pipeline && params.lreads_type == 'pacbio') {
    // loading subreads.bam

    bamfiles = Channel.fromPath(params.pacbio_bamPath).map { file -> tuple(file, file + ".pbi") }
    pacbio_fastq = Channel.empty()

    // loading pacbio fastqs

  } else if (params.longReads && params.run_longreads_pipeline && params.lreads_type == 'pacbio') {
    pacbio_fastq = Channel.fromPath(params.longReads)
    bamfiles = Channel.empty()

    // not loading

  } else {
    bamfiles = Channel.empty()
    pacbio_fastq = Channel.empty()
}

/*

      Loading Pacbio Legacy .bas.h5 data

*/
if (params.lreads_type == 'pacbio' && params.pacbio_h5Path && params.run_longreads_pipeline) {
  String file = params.pacbio_h5Path;
  def h5_name = file.minus(".bas.h5");
  def h5_bax = h5_name + "*.bax.h5"
  println h5_bax ;
  h5bax = Channel.fromPath(h5_bax).collect() ;
  h5files = Channel.fromPath(params.pacbio_h5Path)
} else {
  h5bax = Channel.empty()
  h5files = Channel.empty()
}

/*

          PIPELINE BEGINS

*/

/*
        First Block - Short Reads Pre Processing

        Running with Paired End Short Reads.

*/

/*

        STEP 1 - FastQC

*/

process fastqcPaired {
  publishDir outdir, mode: 'copy',
  // This line saves all the zip files in a folder named "zips"
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
  // This line loads the docker container needed
  container 'fmalmeida/compgen:PREPROCESS'
  tag "Executing FastQC with paired end reads."

    input:
    set val(id), file(sread1), file(sread2) from fastqc_paired

    output:
    file "fastqc_${id}/*_fastqc.{zip,html}"

    when:
    // Execute this process only when having paired end short reads.
    params.reads_size == 2 && params.run_shortreads_pipeline

    script:
    """
    mkdir fastqc_${id} ;
    fastqc -t ${params.threads} -o fastqc_${id} -q $sread1 $sread2
    """
}

/*

        STEP 2 - Trim Galore!

*/

process trimGalorePaired {
    publishDir outdir, mode: 'copy',
        saveAs: {filename ->
    // This line saves the files with specific sufixes in specific folders
            if (filename.indexOf("_fastqc") > 0) "trimGaloreFastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "trimGaloreFastQC/$filename"
            else if (filename.indexOf(".fq.gz") > 0) "zips/TrimGaloreFastQ/$filename"
            else null
        }
    // This line loads the docker container needed
    container 'fmalmeida/compgen:PREPROCESS'
    tag "Executing TrimGalore with paired end reads."

    input:
    set val(id), file(sread1), file(sread2) from trimgalore_paired

    output:
    set val(id), file("${id}_1.fq.gz"), file("${id}_2.fq.gz") into galore_trimmed_paired_reads
    file "*trimming_report.txt"
    file "*_fastqc.{zip,html}"

    when:
    // Execute this process only when having paired end short reads.
    params.reads_size == 2 && params.run_shortreads_pipeline

    script:
    // Loads Optional Parameters
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    TrimQuality = params.quality_trim != null ? params.quality_trim : ''

    """
    trim_galore --paired -q ${params.quality_trim} \
    --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $sread1 $sread2 ;
    mv *_val_1.fq.gz ${id}_1.fq.gz ;
    mv *_val_2.fq.gz ${id}_2.fq.gz ;
    """
}

/*

        STEP 3 - Read Correction

*/

// Checks if lighter is meant to be executed
lighter_input = (params.lighter_execute) ? Channel.empty().mix(galore_trimmed_paired_reads)
                                         : Channel.empty()

process lighterCorrectionPaired {
   publishDir outdir, mode: 'copy',
   saveAs: {filename ->
   // This line saves the files with specific sufixes in specific folders
            if (filename.indexOf(".cor{.fq.gz, .fq}") > 0) "lighterCorrectedOut/$filename"
            else "lighterOutput/$filename"}
   // This line loads the docker container needed
   container 'fmalmeida/compgen:PREPROCESS'
   tag "Executing Ligther (read correction) with paired end reads."

   input:
   set val(id), file(sread1), file(sread2) from lighter_input

   output:
   set val(id), file('*_1.cor.fq.gz'), file('*_2.cor.fq.gz') into lighter_corrected_paired_reads
   file 'fastqc'

   when:
   // Execute this process only when having paired end short reads.
   params.reads_size == 2 && params.run_shortreads_pipeline && params.lighter_execute

   script:
   // Executes the Command with the user's defined alpha parameter
   if (params.lighter_alpha)
   """
   lighter -r $sread1 -r $sread2 -k ${params.lighter_kmer} ${params.lighter_genomeSize} ${params.lighter_alpha} ;
   mkdir fastqc ;
   fastqc -t ${params.threads} -o fastqc -q *_1.cor.fq.gz *_2.cor.fq.gz
   """
   // Lets Lighter automatically find the best alpha parameter
   else
   """
   lighter -r $sread1 -r $sread2 -K ${params.lighter_kmer} ${params.lighter_genomeSize};
   mkdir fastqc ;
   fastqc -t ${params.threads} -o fastqc -q *_1.cor.fq.gz *_2.cor.fq.gz
   """
}

/*

          OPTIONAL STEP - FLASH READ MERGE

*/

// Checks if lighter was executed
flash_input = (params.lighter_execute) ? Channel.empty().mix(lighter_corrected_paired_reads)
                                      : Channel.empty().mix(galore_trimmed_paired_reads)

process flashMerger {
  publishDir outdir, mode: 'copy',
       saveAs: {filename ->
  // This line saves the files with specific sufixes in specific folders
         if (filename.indexOf(".fastq") > 0) "flash_output/$filename"
         else "flash_output/$filename" }
  // This line loads the docker container needed
  container 'fmalmeida/compgen:PREPROCESS'
  tag "Executing FLASH read merger with paired end reads."

  input:
  set val(id), file(sread1), file(sread2) from flash_input

  output:
  file "flash_out*"

  when:
  // Execute this process only when desired and with paired end short reads.
  params.flash_execute && params.run_shortreads_pipeline

  script:
  """
  source activate flash ;
  flash -q -o flash_out -z -t ${params.threads} $sread1 $sread2 &> flash.log;
  """

}

/*

        Running with unpaired short reads.

*/

/*

        STEP 1 - FastQC

*/

process fastqcSingle {
  publishDir outdir, mode: 'copy',
  // This line saves all the zip files in a folder named "zips"
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
  // This line loads the docker container needed
  container 'fmalmeida/compgen:PREPROCESS'
  tag "Executing FastQC with single end reads."

  input:
  file sread  from fastqc_single

  output:
  file "fastqc_${id}/*_fastqc.{zip,html}" mode flatten

  when:
  // Execute this process only when having single end short reads.
  params.reads_size == 1 && params.run_shortreads_pipeline

  script:
  id = sread.getBaseName()
  """
  mkdir fastqc_${id} ;
  fastqc -t ${params.threads} -o fastqc_${id} -q $sread
  """
}

/*

        STEP 2 - Trim Galore!

*/

process trimGaloreSingle {
  publishDir outdir, mode: 'copy',
  saveAs: {filename ->
  // This line saves the files with specific sufixes in specific folders
            if (filename.indexOf("_fastqc") > 0) "trimGaloreFastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "trimGaloreFastQC/$filename"
            else if (filename.indexOf(".fq.gz") > 0) "zips/TrimGaloreFastQ/$filename"
            else null }
  // This line loads the docker container needed
  container 'fmalmeida/compgen:PREPROCESS'
  tag "Executing TrimGalore with single end reads."

  input:
  file sread  from trimgalore_single

  output:
  file('*.{fq.gz,fq}') into galore_trimmed_single_reads mode flatten
  file "*trimming_report.txt"
  file "*_fastqc.{zip,html}"

  when:
  // Execute this process only when having single end short reads.
  params.reads_size == 1 && params.run_shortreads_pipeline

  script:
  // Loads Optional Parameters
  c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
  tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
  """
  trim_galore --fastqc -q ${params.quality_trim} --gzip $c_r1 $tpc_r1 $sread
  """
}

/*

        STEP 3 - Read Correction

*/

process lighterCorrectionSingle {
  publishDir outdir, mode: 'copy',
  saveAs: {filename ->
  // This line saves the files with specific sufixes in specific folders
            if (filename.indexOf(".cor{.fq.gz, .fq}") > 0) "lighterCorrectedOut/$filename"
            else "lighterOutput/$filename"}
  // This line loads the docker container needed
  container 'fmalmeida/compgen:PREPROCESS'
  tag "Executing Lighter (read correction) with single end reads."

  input:
  file sread  from galore_trimmed_single_reads

  output:
  file "*"

  when:
  // Execute this process only when having single end short reads.
  params.reads_size == 1 && params.run_shortreads_pipeline && params.lighter_execute

  script:
  // Executes the Command with the user's defined alpha parameter
  if (params.lighter_alpha)
  """
  lighter -r $sread -k ${params.lighter_kmer} ${params.lighter_genomeSize} ${params.lighter_alpha};
  mkdir fastqc ;
  fastqc -t ${params.threads} -o fastqc -q *.cor{.fq.gz, .fq}
  """
  // Lets Lighter automatically find the best alpha parameter
  else
  """
  lighter -r $sread -K ${params.lighter_kmer} ${params.lighter_genomeSize};
  mkdir fastqc ;
  fastqc -t ${params.threads} -o fastqc -q *.cor{.fq.gz, .fq}
  """
}

/*
 *      Second Block - Long reads pre-processing
 */

/*
 *      STEP 1 - Extract fasta from h5 files - PACBIO ONLY
 *      Only if user gives h5 data.
 */

process subreadsBamToFastq {
  publishDir outdir, mode: 'copy'
  // Loads the necessary Docker image
  container 'fmalmeida/compgen:PREPROCESS'
  tag "Extracting FASTQ from pacbio .subreads.bam files"

  input:
  set file(input), file(pbi) from bamfiles

  output:
  file "*.fastq" into pacbio_trimmed mode flatten

  when:
  // Sets execution condition. Only when user have pacbio bam data.
  params.lreads_type == 'pacbio' && params.pacbio_bamPath && params.run_longreads_pipeline

  script:
  if (params.lreads_is_barcoded)
  """
  source activate pbtools ;
  bam2fastq -o ${input.baseName} -u --split-barcodes $input ;
  """
  else
  """
  source activate pbtools ;
  bam2fastq -o ${input.baseName} -u $input
  """
}

process legacyH5ToFastq {
  publishDir outdir, mode: 'copy'
  // Loads the necessary Docker image
  container 'fmalmeida/compgen:PREPROCESS'
  tag "Extracting FASTQ from h5 files"

  input:
  file(inputh5) from h5files
  file h5bax

  output:
  file "${inputh5.baseName}.fastq" into pacbio_legacy_extracted mode flatten

  when:
  // Sets execution condition. Only when user have pacbio h5 data.
  params.lreads_type == 'pacbio' && params.pacbio_h5Path && params.run_longreads_pipeline

  script:
  """
  bash5tools.py --outFilePrefix ${inputh5.baseName} --readType subreads \
  --outType fastq --minLength 200 ${inputh5} ;
  """
}

/*

        STEP 2 - Trimming long reads

*/

process porechopTrimming {
  publishDir outdir, mode: 'copy'
  // Loads the necessary Docker image
  container 'fmalmeida/compgen:PREPROCESS'
  tag "Trimming with Porechop"

  input:
  file reads from lreads

  output:
  file "${reads.baseName}_trimmed.fastq" into ont_normal_trimmed mode flatten optional true
  file "porechop_barcodes/*.fastq" into ont_barcodes_trimmed mode flatten optional true

  when:
  // Sets execution condition.
  params.lreads_type == 'nanopore' && params.run_longreads_pipeline

  script:
  if (params.lreads_is_barcoded)
  """
  porechop -i ${reads} -b porechop_barcodes --barcode_threshold 85
  """
  else
  """
  porechop -i ${reads} -t ${params.threads} --format fastq -o ${reads.baseName}_trimmed.fastq ;
  """
}

// Create nanopack input channel

ont_trimmed_reads = Channel.empty().mix(ont_normal_trimmed, ont_barcodes_trimmed)
nanopack_input = Channel.empty().mix(pacbio_fastq, ont_trimmed_reads, pacbio_trimmed, pacbio_legacy_extracted)

/*

        STEP 3 - Getting long reads Sequencing quality

*/

process nanopack {
  publishDir outdir, mode: 'copy'
  // Loads the necessary Docker image
  container 'fmalmeida/compgen:PREPROCESS'
  tag "Checking ONT reads stats with NanoPack"

  input:
  file input from nanopack_input

  output:
  file "*"

  when:
  // Sets execution condition. Only when user wants porechop trimming
  params.run_longreads_pipeline

  script:
  """
  source activate nanopack;
  # Plotting
  NanoPlot -t 3 --fastq $input -o ${input.baseName}_nanoplot -f svg --N50 \
  --title "Kp31 sample" --plots hex dot pauvre kde ;

  # Checking Quality
  nanoQC -o ${input.baseName}_nanoQC $input ;

  # Generate Statistics Summary
  NanoStat --fastq $input -n ${input.baseName}.txt --outdir ${input.baseName}_stats ;
  """
}

// Completition message

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
}

// Header log info
log.info "============================================"
log.info " Docker-based reads pre-processing Pipeline "
log.info "============================================"
def summary = [:]
summary['Output dir']   = params.outDir
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="
