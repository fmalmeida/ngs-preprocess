# ngs-preprocess pipeline changelog

The tracking for changes started in v2.2

## v2.7.1

* [[#36]](https://github.com/fmalmeida/ngs-preprocess/issues/36) - Include nf-tests to the pipeline.
* [[#40]](https://github.com/fmalmeida/ngs-preprocess/issues/40) - Fix a problem when downloading PACBIO data.
* Adjust short-reads workflow filter to allow analysis from SRA for platform key-values: illumina, bgiseq and dnbseq

## v2.7.0 -- [2024-Apr-30]

* [[#34]](https://github.com/fmalmeida/ngs-preprocess/issues/34) - Included a new tool in the pipeline (`porechop ABI`). Removed `tracedir` parameter. Build new docker image.

## v2.6.3 -- [2024-Fev-23]

* [[#32](https://github.com/fmalmeida/ngs-preprocess/issues/32)] - Add as output a template samplesheet that can be readily used as input for MpGAP to assemble each downloaded read on its own.

## v2.6.2 -- [2024-Jan-19]

* [[#24](https://github.com/fmalmeida/ngs-preprocess/issues/24)] - Added documentation on generated outputs, as requested in paper review
* [[#25](https://github.com/fmalmeida/ngs-preprocess/issues/25)] - Added documentation exemplifying use for non bacterial dataset, as requested in paper review
* Small fixes on SRA module to allow downloading bgiseq data, and avoid running on empty lines
* Updated Dockerfile to name the maximum of versions possible

> Nothing has changed in terms of how tools are called and used, thus the docker image still the same. In fact, patch/fix releases (x.x.x) will always use the docker from breaking/features release (x.x)

## v2.6.1 -- [2023-Oct-26]

* [[#21](https://github.com/fmalmeida/ngs-preprocess/issues/21)] - Added/updated citation information
* Small fix on pacbio download from sra

> Nothing has changed in terms of how tools are called and used, thus the docker image still the same. In fact, patch/fix releases (x.x.x) will always use the docker from breaking/features release (x.x)

## v2.6 -- [2023-Apr-16]

Small bug fixes on SRA_FETCH download modules. Updated pipeline to understand Illumina / BGIseq sequences for short reads. Updated Docker image with named versions of tools whenever possible.

> More tools have been added so the versioning and docker image have now changed to v2.6.

## v2.5 -- [2022-Oct-30]

Add possibility for users to automatically fetch fastq files from SRA NCBI database. For that, users just need to use the `--sra_ids` parameter, passing a file with a list of SRA RunIDs, one per line.

> More tools have been added so the versioning and docker image have now changed to v2.5.

## v2.4.2 -- [2022-Oct-17]

Cleanup change. Short reads output are are now written as "preprocessed_reads/short_reads" instead of "preprocessed/illumina" as sometimes other technology may be used.

> Nothing has changed in terms of how tools are called and used, thus the docker image still the same. In fact, patch/fix releases (x.x.x) will always use the docker from breaking/features release (x.x)

## v2.4.1 -- [2022-Feb-21]

This version addresses the changes discussed in [issue #10](https://github.com/fmalmeida/ngs-preprocess/issues/10). It has three main changes:

1. Added standard NF allocation resource rules as it is done by nf-core community
    * It also uses templates of CLI help and logging messages from nf-core community.
2. Re-organized config files to keep structure cleaner
3. Changed the standar profile which will not load docker by default anymore. As it is the common practice for NF pipelines, user must explicitily select between docker/conda/singularity profiles.

> Nothing has changed in terms of how tools are called and used, thus the docker image still the same. In fact, patch/fix releases (x.x.x) will always use the docker from breaking/features release (x.x)

## v2.4

This version addresses the changes discussed in [issue #7](https://github.com/fmalmeida/ngs-preprocess/issues/7). It has three main changes:

1. The directory of outputs have been reorganized and the output files extension have been standardized in order to make it easier and simpler to use this pipeline as the first step of many analyses. It must look like this now:

```console
{OUTPUT}
├── final_output  # the final preprocessed reads of any tech (organized in folders)
│   ├── illumina  # final preprocessed illumina reads
│   │   ├── SRR8482585_30X.merged.fq.gz
│   │   ├── SRR8482585_30X_R1.unmerged.fq.gz
│   │   └── SRR8482585_30X_R2.unmerged.fq.gz
│   ├── nanopore  # final preprocessed nanopore reads
│   │   └── kp_30X.filtered.fq.gz
│   └── pacbio    # final preprocessed pacbio reads
│       └── lima.bc1018--bc1018.filtered.fq.gz
└── preprocessing_outputs   # here will be saved (per tech) the files (QC, Logs, Stats) generated during preprocessing
    ├── illumina
    │   ├── SRR8482585_30X_fastp.html
    │   └── SRR8482585_30X_fastp.json
    ├── nanopore
    │   ├── QC
    │   │   ├── kp_30X_nanoQC
    │   │   ├── kp_30X_nanoplot
    │   │   └── kp_30X_stats
    │   └── porechop
    │       └── kp_30X.trimmed.fq.gz
    └── pacbio
        ├── QC
        │   ├── lima.bc1018--bc1018_nanoQC
        │   ├── lima.bc1018--bc1018_nanoplot
        │   └── lima.bc1018--bc1018_stats
        └── bam2fastq
            ├── lima.bc1018--bc1018.bam.pbi
            └── lima.bc1018--bc1018.fq.gz
```
2. Secondly, the pipelines tools/dependencies for preprocessing Illumina reads have changed. Now, instead of using flash, trim_galore, fastqc and lighter, the pipeline will use only [fastp](https://github.com/OpenGene/fastp) which performs all the steps that were done by these tools, but much faster since it is written in C++.
3. The parameter `--outdir` is now `--output` for better readability

## v2.3

This version is related to issue https://github.com/fmalmeida/ngs-preprocess/issues/5. It re-configured the pipeline to be more likely to nf-core pipelines, enabling users to run it using docker, conda or singularity. More explanation on how to run with different profiles is given at: https://github.com/fmalmeida/ngs-preprocess/tree/master#selecting-between-profiles

Moreover, some small fixes have been done in all the workflows so that the pre-processed reads are stored in sub-folders named as the reads basename. Before the pipeline was placing it all under the `--outdir` directory without sub-folders organization.

## v2.2

This version have a few additions to the pipeline workflow, they are highlighted and explained below:

### nf-core schema

We have added a nextflow parameter schema in json that is compliant with nf-core. This enables that users trigger the graphical interface for configuration and execution of the pipeline via [nf-core launch](https://nf-co.re/launch) utility.

```bash
# It is triggered as
nf-core launch fmalmeida/ngs-preprocess
```

### pacbio hifi

We have added the possibility for users to call high fidelity reads with [pacbio ccs](https://ccs.how/) software. For this, users must use the parameter `--pacbio_get_hifi`

### pacbio demultiplexing

We have changed the demultiplexing process to use [lima](https://github.com/PacificBiosciences/barcoding) instead of [bam2fastx](https://github.com/PacificBiosciences/bam2fastx). In order to demultiplex the reads, users must give the pacbio barcodes (fasta or xml) via `--pacbio_barcodes` parameter.

Be careful, lima enables that users demultiplex data by same, different or any kind of association of barcodes. By default the pipeline checks for reads that have the same barcodes. Checkout the `--pacbio_barcode_design` parameter.

### longreads filtering

Using [nanofilt](https://github.com/wdecoster/nanofilt) we have added the possibility that users filter their long reads (ONT or Pacbio) based on min. length and quality thresholds. Check the `--lreads_min_length` and `--lreads_min_quality` parameters.
