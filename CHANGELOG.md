# ngs-preprocess pipeline changelog

The tracking for changes started in v2.2

## v2.2

This versions have a few additions to the pipeline workflow, they are highlighted and explained below:

### nf-core schema

We have added a nextflow parameter schema in json that is compliant with nf-core. This enables that users trigger the graphical interface for configuration and execution of the pipeline via [nf-core launch](https://nf-co.re/launch) utility, also it is possible to track the pipeline execution with [nextflow tower](https://tower.nf/).

```bash
# It is triggered as
nf-core launch fmalmeida/ngs-preprocess
```

Checkout the paremeters `--use_tower` and `--tower_token` to activate pipeline execution in nextflow tower.

### pacbio hifi

We have added the possibility for users to call high fidelity reads with [pacbio ccs](https://ccs.how/) software. For this, users must use the parameter `--pacbio_get_hifi`

### pacbio demultiplexing

We have changed the demultiplexing process to use [lima](https://github.com/PacificBiosciences/barcoding) instead of [bam2fastx](https://github.com/PacificBiosciences/bam2fastx). In order to demultiplex the reads, users must give the pacbio barcodes (fasta or xml) via `--pacbio_barcodes` parameter.

Be careful, lima enables that users demultiplex data by same, different or any kind of association of barcodes. By default the pipeline checks for reads that have the same barcodes. Checkout the `--pacbio_barcode_design` parameter.

### longreads filtering

Using [nanofilt](https://github.com/wdecoster/nanofilt) we have added the possibility that users filter their long reads (ONT or Pacbio) based on min. length and quality thresholds. Check the `--lreads_min_length` and `--lreads_min_quality` parameters.
