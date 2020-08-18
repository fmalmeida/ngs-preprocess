# ngs-preprocess pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3634044.svg)](https://doi.org/10.5281/zenodo.3634044) ![](https://img.shields.io/github/v/release/fmalmeida/ngs-preprocess) [![Build Status](https://travis-ci.com/fmalmeida/ngs-preprocess.svg?branch=master)](https://travis-ci.com/fmalmeida/ngs-preprocess) ![](https://img.shields.io/docker/cloud/build/fmalmeida/ngs-preprocess) [![Documentation Status](https://readthedocs.org/projects/ngs-preprocess/badge/?version=latest)](https://ngs-preprocess.readthedocs.io/en/latest/?badge=latest) ![](https://img.shields.io/badge/Nextflow-v20.07-yellowgreen)


ngs-preprocess is an easy to use nextflow docker-based pipeline that uses state-of-the-art software for pre-procesing ngs reads of Illumina, Pacbio and Oxford Nanopore Technologies and has only two dependencies: [Docker](https://www.docker.com/) and [Nextflow](https://github.com/nextflow-io/nextflow). It wraps up the following software:

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
* [FLASH](https://ccb.jhu.edu/software/FLASH/)
* [Lighter](https://github.com/mourisl/Lighter)
* [Porechop](https://github.com/rrwick/Porechop)
* [pycoQC](https://github.com/a-slide/pycoQC)
* [bam2fastx](https://github.com/PacificBiosciences/bam2fastx)
* [bax2bam](https://github.com/PacificBiosciences/bax2bam)
* [NanoPack](https://github.com/wdecoster/nanopack).

## Table of contents

* [Requirements](https://github.com/fmalmeida/ngs-preprocess#requirements)
* [Quickstart](https://github.com/fmalmeida/ngs-preprocess#quickstart)
* [Documentation](https://github.com/fmalmeida/ngs-preprocess#documentation)
  * [Full usage](https://github.com/fmalmeida/ngs-preprocess#usage)
  * [Usage Examples](https://github.com/fmalmeida/ngs-preprocess#usage-examples)
  * [Configuration File](https://github.com/fmalmeida/ngs-preprocess#using-the-configuration-file)

## Requirements

* Unix-like operating system (Linux, macOS, etc)
* Nextflow (version 20.01 or higher)
* Java 8
* Docker
  * Image: `fmalmeida/ngs-preprocess`

## Quickstart

1. If you don't have it already install Docker in your computer. Read more [here](https://docs.docker.com/).
    * You can give this [in-house script](https://github.com/fmalmeida/bioinfo/blob/master/dockerfiles/docker_install.sh) a try.
    * After installed, you need to download the required Docker image

          docker pull fmalmeida/ngs-preprocess

2. Install Nextflow (version 20.01 or higher):

       curl -s https://get.nextflow.io | bash

3. Give it a try:

       nextflow fmalmeida/ngs-preprocess --help

## Documentation

### Usage

Checkout the full usage help with `nextflow run fmalmeida/ngs-preprocess --help`. Users are encouraged to read the [complete documentation](https://ngs-preprocess.readthedocs.io/en/latest/?badge=latest).

### Usage examples:

> Illumina paired end reads. Since it will always be a pattern match, example "illumina/SRR9847694_{1,2}.fastq.gz", it MUST ALWAYS be double quoted as the example below.

    ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir outputs/illumina_paired \
    --shortreads "illumina/SRR9847694_{1,2}.fastq.gz" --shortreads_type "paired" \
    --lighter_execute --lighter_genomeSize 4.6m --clip_r1 5 --three_prime_clip_r1 5 \
    --clip_r2 5 --three_prime_clip_r2 5 --quality_trim 30 --flash_execute

> Illumina single end reads. Multiple files at once, using fixed number of bases to be trimmed. If multiple unpaired reads are given as input at once, pattern MUST be double quoted: "SRR9696*.fastq.gz"

    ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir outputs/illumina_single \
    --shortreads "illumina/SRR9696*.fastq.gz" --shortreads_type "single" \
    --clip_r1 5 --three_prime_clip_r1 5

> ONT barcoded reads (fastq)

    ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir outputs/ont \
    --nanopore_fastq kpneumoniae_25X.fastq --nanopore_is_barcoded

> Pacbio raw barcoded (subreads.bam) reads

    ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outdir outputs/pacbio \
    --pacbio_bamPath pacbio/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreads.bam \
    --pacbio_is_barcoded

> Pacbio raw (legacy .bas.h5 to subreads.bam) reads

    ./nextflow run fmalmeida/ngs-preprocess --pacbio_h5Path E01_1/Analysis_Results/ \
    --outdir E01_1/Analysis_Results/preprocessed --threads 3

> Running with a configuration file

    ./nextflow run fmalmeida/ngs-preprocess -c nextflow.config

## Using the configuration file

All the parameters showed above can be, and are advised to be, set through the configuration file. When a configuration file is set the pipeline is run by simply executing `nextflow run fmalmeida/ngs-preprocess -c ./configuration-file`

Your configuration file is what will tell to the pipeline the type of data you have, and which processes to execute. Therefore, it needs to be correctly set up.

Create a configuration file in your working directory:

* Complete config:

      nextflow run fmalmeida/ngs-preprocess --get_full_config

* For Illumina data:

      nextflow run fmalmeida/ngs-preprocess --get_illumina_config

* For Pacbio data:

      nextflow run fmalmeida/ngs-preprocess --get_pacbio_config

* For ONT data:

      nextflow run fmalmeida/ngs-preprocess --get_ont_config

# Citation

To cite this tool please refer to our Zenodo tag or directly via the github url.

Users are encouraged to cite the programs used in this pipeline whenever they are used. They are: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [TrimGalore](https://github.com/FelixKrueger/TrimGalore), [FLASH](https://ccb.jhu.edu/software/FLASH/), [Lighter](https://github.com/mourisl/Lighter), [Porechop](https://github.com/rrwick/Porechop), [pycoQC](https://github.com/a-slide/pycoQC), [bax2bam](https://github.com/PacificBiosciences/bax2bam), [bam2fastq](https://github.com/PacificBiosciences/bam2fastx) and [NanoPack](https://github.com/wdecoster/nanopack).
