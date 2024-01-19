# Output files

Here, using the results produced in the [Non-bacterial dataset section](non_bacteria.md#), we give users a glimpse over the main outputs produced by bacannot. The command used in the quickstart wrote the results under the `preprocessed_reads` directory.

!!! note

    Please take note that the pipeline uses the directory set with the `--output` parameter as a storage place in which it will create a folder for the final pre-processed reads and for the intermediate files, separated by sequencing technology.

## Directory tree

After a successful execution, you will have something like this:

```bash

# Directory tree from the running dir
preprocessed_reads
├── final_output                                              # directory containing the final results of the data cleaning                   
│   └── nanopore
│       └── SRR23337893.filtered.fq.gz
├── pipeline_info                                             # directory containing the nextflow execution reports
│   ├── ngs_preprocess_report_2023-11-18_10-07-36.html
│   ├── ngs_preprocess_timeline_2023-11-18_10-07-36.html
│   ├── ngs_preprocess_tracing_2023-11-18_10-07-36.txt
├── preprocessing_outputs                                     # directory containing the intermediate files produced by the tools used during pre-processing, and, QC
│   └── nanopore
│       ├── porechop
│       └── QC
└── SRA_FETCH                                                 # directory containing the intermediate files when downloading data from SRA
    ├── FASTQ
    │   └── SRR23337893_data
    └── SRR23337893_sra_runInfo.csv
```
