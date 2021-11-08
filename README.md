<img src="images/lOGO_3.png" width="300px">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3451405.svg)](https://doi.org/10.5281/zenodo.3451405)
[![Releases](https://img.shields.io/github/v/release/fmalmeida/ngs-preprocess)](https://github.com/fmalmeida/ngs-preprocess/releases)
[![Documentation](https://img.shields.io/badge/Documentation-readthedocs-brightgreen)](https://ngs-preprocess.readthedocs.io/en/latest/?badge=latest)
[![Dockerhub](https://img.shields.io/badge/Docker-fmalmeida/ngs--preprocess-informational)](https://hub.docker.com/r/fmalmeida/ngs-preprocess)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/fmarquesalmeida)
[![License](https://img.shields.io/badge/License-GPL%203-black)](https://github.com/fmalmeida/ngs-preprocess/blob/master/LICENSE)

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

ngs-preprocess is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. It is an easy to use pipeline that uses state-of-the-art software for quality check and pre-processing ngs reads of Illumina, Pacbio and Oxford Nanopore Technologies.

It wraps up the following software:

| Step | tools |
| :--- | :---- |
| Illumina pre-processing | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [TrimGalore](https://github.com/FelixKrueger/TrimGalore), [FLASH](https://ccb.jhu.edu/software/FLASH/), [Lighter](https://github.com/mourisl/Lighter) |
| Nanopore pre-processing | [Porechop](https://github.com/rrwick/Porechop), [pycoQC](https://github.com/tleonardi/pycoQC), [NanoPack](https://github.com/wdecoster/nanopack) |
| Pacbio pre-processing | [bam2fastx](https://github.com/PacificBiosciences/bam2fastx), [bax2bam](https://github.com/PacificBiosciences/bax2bam), [lima](https://github.com/PacificBiosciences/barcoding), [pacbio ccs](https://ccs.how/) |

## Further reading

This pipeline has two complementary pipelines (also written in nextflow) for [genome assembly](https://github.com/fmalmeida/mpgap) and [prokaryotic genome annotation](https://github.com/fmalmeida/bacannot) that can give the user a complete workflow for bacterial genomics analyses.

## Quickstart

1. Install Nextflow (version 20.01 or higher):
    + `curl -s https://get.nextflow.io | bash`
2. Give it a try:
    + `nextflow run fmalmeida/ngs-preprocess --help`
3. Start running your analysis
    
    ```console
    nextflow run fmalmeida/ngs-preprocess -profile <docker/singularity/conda>
    ```

:fire: Please read the documentation below on [selecting between conda, docker or singularity](https://github.com/fmalmeida/ngs-preprocess/tree/master#selecting-between-conda-docker-and-singularity) profiles, since the tools will be made available differently depending on the profile desired.

## Documentation

### Selecting between profiles

Nextflow profiles are a set of "sensible defaults" for the resource requirements of each of the steps in the workflow, that can be enabled with the command line flag `-profile`. You can learn more about nextflow profiles at:

+ https://nf-co.re/usage/configuration#basic-configuration-profiles
+ https://www.nextflow.io/docs/latest/config.html#config-profiles

The pipeline have "standard profiles" set to run the workflows with either conda, docker or singularity using the [local executor](https://www.nextflow.io/docs/latest/executor.html), which is nextflow's default and basically runs the pipeline processes in the computer where Nextflow is launched. If you need to run the pipeline using another executor such as sge, lsf, slurm, etc. you can take a look at [nextflow's manual page](https://www.nextflow.io/docs/latest/executor.html) to proper configure one in a new custom profile set in your personal copy of [ngs-preprocess config file](https://github.com/fmalmeida/ngs-preprocess/blob/master/nextflow.config) and take advatage that nextflow allows multiple profiles to be used at once, e.g. `-profile conda,sge`.

By default, if no profile is chosen, the pipeline will "load the docker profile". Available pre-set profiles for this pipeline are: docker, conda, singularity, you can choose between them as follows:

* conda

    ```console
    nextflow run fmalmeida/ngs-preprocess -profile conda [options]
    ```

* docker
    
    ```console
    nextflow run fmalmeida/ngs-preprocess -profile docker [options]
    ```

* singularity
    
    ```console
    nextflow run fmalmeida/ngs-preprocess -profile singularity [options]
    ```

:note: Please use conda as last resource since the packages will not be packed, problems may arise, and nextflow will trigger an installation every time which may consume plenty of time.

### Usage

<a href="https://ngs-preprocess.readthedocs.io/en/latest/index.html"><strong>Users are advised to read the complete documentation »</strong></a>

* Complete command line explanation of parameters:
    + `nextflow run fmalmeida/ngs-preprocess --help`
* See usage examples in the command line:
    + `nextflow run fmalmeida/ngs-preprocess --examples`

### Command line usage examples

Command line executions are exemplified [in the manual](https://ngs-preprocess.readthedocs.io/en/latest/examples.html).

**Remember**: Whenever using REGEX for a pattern match, for example "illumina/SRR9847694_{1,2}.fastq.gz", it MUST ALWAYS be inside double quotes.

### Using the configuration file

All the parameters showed above can be, and are advised to be, set through the configuration file. When a configuration file is set the pipeline is run by simply executing `nextflow run fmalmeida/ngs-preprocess -c ./configuration-file`

Your configuration file is what will tell to the pipeline the type of data you have, and which processes to execute. Therefore, it needs to be correctly set up.

Create a configuration file in your working directory:

```console
nextflow run fmalmeida/ngs-preprocess [ --get_full_config ] [ --get_illumina_config ] [ --get_pacbio_config ] [ --get_ont_config ]
```

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
