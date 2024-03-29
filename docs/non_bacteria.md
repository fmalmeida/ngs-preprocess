# Example use case of a non-bacterial dataset

As part of the reviews received in the [published paper](https://doi.org/10.12688/f1000research.139488.1) of this pipeline, it was requested to provide an example use case of the preprocessing (this pipeline) and the assembly pipeline ([MpGAP](https://github.com/fmalmeida/MpGAP)) using a non-bacterial dataset in order to demonstrate that these two pipelines are **not** specific to bacterial genomes.

Thus, this is what we are providing in this page. The preprocessed dataset produced here, will be the data used in the assembly pipeline as a way of showing it connection.

## Get the data

We are going to analyse the dataset `SRR23337893` of an _Aspergillus fumigatus_ isolate. The pipeline is capable of automatically downloading the data from SRA, thus, we simply need to generate an input file that contains a list of SRA IDs to be analysed, in our case, only this one.


```bash
echo 'SRR23337893' >  aspergillus_fumigatus_sra.txt
```

This is the only input required to run the pipeline

## Preprocessing the data

Outputs will be at `preprocessed_reads`.

```bash
nextflow run fmalmeida/ngs-preprocess \
    -profile docker \
    --sra_ids aspergillus_fumigatus_sra.txt \
    --output ./preprocessed_reads \
    --max_cpus 10 \
    --max_memory 20.GB \
    --lreads_min_length 750 \
    --lreads_min_quality 10
```

## Afterwards

The generated outputs are displayed as an example of the outputs generated by the pipeline in the [outputs](outputs.md#) page.

Finally, the preprocessed reads generated in this example, are the ones that are used in the "non-bacterial" example use case displayed in the assembly pipeline, [MpGAP](https://mpgap.readthedocs.io/en/latest/index.html).
