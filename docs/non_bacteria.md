# Example use case of a non-bacterial dataset

As part of the reviews received in the [published paper](https://doi.org/10.12688/f1000research.139488.1) of this pipeline, it was requested to provide an example use case of the preprocessing (this pipeline) and the assembly pipeline ([MpGAP](https://github.com/fmalmeida/MpGAP)) using a non-bacterial dataset in order to demonstrate that these two pipelines are **not** specific to bacterial genomes.

Thus, this is what we are providing in this page. The preprocessed dataset produced here, will be the data used in the assembly pipeline as a way of showing it connection.

## Get the data

We have made this subsampled dataset available in [Figshare](https://figshare.com/articles/dataset/Illumina_pacbio_and_ont_sequencing_reads/14036585).


```bash
# Download data from figshare
wget -O reads.zip https://ndownloader.figshare.com/articles/14036585/versions/4

# Unzip
unzip reads.zip
```

Now we have the necessary data to perform the quickstart.

!!! note "Where my outputs go?"

    The pipeline will always use the fastq file name as prefix for sub-folders and output files. For instance, if users use a fastq file named SRR7128258.fastq the output files and directories will have the string "SRR7128258" in it.

## Preprocessing the data

Outputs will be at `preprocessed_reads`.

!!! warning "Watch your REGEX"

    Whenever using REGEX for a pattern match, for example "illumina/SRR9847694_{1,2}.fastq.gz" or "illumina/SRR*.fastq.gz", it MUST ALWAYS be inside double quotes.

    **Remember:** the pipeline does not concatenate the reads. Whenever you use a pattern such as \* with unpaired reads the pipeline will process each read separately.

```bash
# Running for both illumina and nanopore data
nextflow run fmalmeida/ngs-preprocess \
    -profile docker \
    --output preprocessed_reads \
    --max_cpus 4 \
    --shortreads "SRR8482585_30X_{1,2}.fastq.gz" \
    --shortreads_type "paired" \
    --fastp_correct_pairs \
    --fastp_merge_pairs \
    --nanopore_fastq "SRX5299443_30X.fastq.gz" \
    --lreads_min_length 1000 \
    --lreads_min_quality 10
```

## Using test profile


As for version v2.5, users can also used a pre-configured test profile which will automatically load a list of SRA run ids for download.

```bash
# Running for both short and long reads data
nextflow run fmalmeida/ngs-preprocess -profile docker,test
```

## Afterwards

Now you can used these datasets to, for example, assemble and annotate a genome. For this, check out the [MpGAP](https://mpgap.readthedocs.io/en/latest/index.html) and [Bacannot](https://bacannot.readthedocs.io/en/latest/index.html) pipelines that we've developed for such tasks.
