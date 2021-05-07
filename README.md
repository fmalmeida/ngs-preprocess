[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3451405.svg)](https://doi.org/10.5281/zenodo.3451405) [![Releases](https://img.shields.io/github/v/release/fmalmeida/ngs-preprocess)](https://github.com/fmalmeida/ngs-preprocess/releases) [![Documentation](https://img.shields.io/badge/Documentation-readthedocs-brightgreen)](https://ngs-preprocess.readthedocs.io/en/latest/?badge=latest) [![Dockerhub](https://img.shields.io/badge/Docker-fmalmeida/ngs--preprocess-informational)](https://hub.docker.com/r/fmalmeida/ngs-preprocess) [![Docker build](https://img.shields.io/docker/cloud/build/fmalmeida/ngs-preprocess)](https://hub.docker.com/r/fmalmeida/ngs-preprocess) ![Docker Pulls](https://img.shields.io/docker/pulls/fmalmeida/ngs-preprocess) [![Nextflow version](https://img.shields.io/badge/Nextflow%20>=-v20.07-important)](https://www.nextflow.io/docs/latest/getstarted.html) [![License](https://img.shields.io/badge/License-GPL%203-black)](https://github.com/fmalmeida/ngs-preprocess/blob/master/LICENSE)

<p align="center">
  <!-- <a href="https://github.com/othneildrew/Best-README-Template">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a> -->

  <h1 align="center">ngs-preprocess pipeline</h2>

  <p align="center">
    <h3 align="center">A pipeline for preprocessing short and long sequencing reads</h3>
    <br />
    <a href="https://ngs-preprocess.readthedocs.io/en/latest/index.html"><strong>See the documentation »</strong></a>
    <br />
    <br />
    <a href="https://github.com/fmalmeida/ngs-preprocess/issues">Report Bug</a>
    ·
    <a href="https://github.com/fmalmeida/ngs-preprocess/issues">Request Feature</a>
  </p>
</p>

## About

ngs-preprocess is an easy to use nextflow docker-based pipeline that uses state-of-the-art software for pre-procesing ngs reads of Illumina, Pacbio and Oxford Nanopore Technologies and has only two dependencies: [Docker](https://www.docker.com/) and [Nextflow](https://github.com/nextflow-io/nextflow). It wraps up the following software:

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
* [FLASH](https://ccb.jhu.edu/software/FLASH/)
* [Lighter](https://github.com/mourisl/Lighter)
* [Porechop](https://github.com/rrwick/Porechop)
* [pycoQC](https://github.com/a-slide/pycoQC)
* [bam2fastx](https://github.com/PacificBiosciences/bam2fastx)
* [bax2bam](https://github.com/PacificBiosciences/bax2bam)
* [lima](https://github.com/PacificBiosciences/barcoding)
* [pacbio ccs](https://ccs.how/)
* [NanoPack](https://github.com/wdecoster/nanopack).

## Further reading

This pipeline has two complementary pipelines (also written in nextflow) for [genome assembly](https://github.com/fmalmeida/mpgap) and [prokaryotic genome annotation](https://github.com/fmalmeida/bacannot) that can give the user a complete workflow for bacterial genomics analyses.

## Requirements

* Unix-like operating system (Linux, macOS, etc)
* Nextflow (version 20.01 or higher)
* Java 8
* Docker
  * Image: `fmalmeida/ngs-preprocess`

## Quickstart

1. If you don't have it already install [Docker](https://docs.docker.com/) in your computer.
    * After installed, you need to download the required Docker image

          docker pull fmalmeida/ngs-preprocess

> Each release is accompanied by a Dockerfile in the docker folder. When using older releases, users can create the correct image using
the Dockerfile that goes alongside with the release (Remember to give the image the correct name, as it is in dockerhub and the nextflow script).
The latest release will always have its docker image in dockerhub.

2. Install Nextflow (version 20.01 or higher):

       curl -s https://get.nextflow.io | bash

3. Give it a try:

       nextflow run fmalmeida/ngs-preprocess --help

> Users can get let the pipeline always updated with: `nextflow pull fmalmeida/ngs-preprocess`

## Documentation

### Usage

<a href="https://ngs-preprocess.readthedocs.io/en/latest/index.html"><strong>Users are advised to read the complete documentation »</strong></a>

* Complete command line explanation of parameters:
    + `nextflow run fmalmeida/ngs-preprocess --help`
* See usage examples in the command line:
    + `nextflow run fmalmeida/ngs-preprocess --examples`

### Command line usage examples

Command line executions are exemplified [in the manual](https://ngs-preprocess.readthedocs.io/en/latest/examples.html).

### Using the configuration file

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

### Interactive graphical configuration and execution

#### Via NF tower launchpad (good for cloud env execution)

Nextflow has an awesome feature called [NF tower](https://tower.nf). It allows that users quickly customise and set-up the execution and configuration of cloud enviroments to execute any nextflow pipeline from nf-core, github (this one included), bitbucket, etc. By having a compliant JSON schema for pipeline configuration it means that the configuration of parameters in NF tower will be easier because the system will render an input form.

Checkout more about this feature at: https://seqera.io/blog/orgs-and-launchpad/

<p align="center">
<img src="https://j.gifs.com/GRnqm7.gif" width="500px"/>
</p>

#### Via nf-core launch (good for local execution)

Users can trigger a graphical and interactive pipeline configuration and execution by using [nf-core launch](https://nf-co.re/launch) utility. nf-core launch will start an interactive form in your web browser or command line so you can configure the pipeline step by step and start the execution of the pipeline in the end.

```bash
# Install nf-core
pip install nf-core

# Launch the pipeline
nf-core launch fmalmeida/ngs-preprocess
```

It will result in the following:

<p align="center">
<img src="./images/nf-core-asking.png" width="500px"/>
</p>

<p align="center">
<img src="./images/nf-core-gui.png" width="400px"/>
</p>

# Citation

To cite this tool please refer to our Zenodo tag or directly via the github url.

Users are encouraged to cite the programs used in this pipeline whenever they are used. They are: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [TrimGalore](https://github.com/FelixKrueger/TrimGalore), [FLASH](https://ccb.jhu.edu/software/FLASH/), [Lighter](https://github.com/mourisl/Lighter), [Porechop](https://github.com/rrwick/Porechop), [pycoQC](https://github.com/a-slide/pycoQC), [bax2bam](https://github.com/PacificBiosciences/bax2bam), [bam2fastq](https://github.com/PacificBiosciences/bam2fastx), [lima](https://github.com/PacificBiosciences/barcoding), [pacbio ccs](https://ccs.how/) and [NanoPack](https://github.com/wdecoster/nanopack).
