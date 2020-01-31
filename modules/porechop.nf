#!/usr/bin/env nextflow

/*
 * Definition of porechop module
 */
process porechop {
  publishDir "/home/falmeida/Downloads/porechopTESTE", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Trimming longreads with Porechop"

  input:
    file reads
    val threads
    val barcode
  output:
    file "${reads.baseName}_trimmed.fastq" optional true
    file "porechop_barcodes/*.fastq" optional true

  script:
    if (barcode)
    """
    porechop -i ${reads} -b porechop_barcodes --barcode_threshold 85
    """
    else
    """
    porechop -i ${reads} -t ${threads} --format fastq -o ${reads.baseName}_trimmed.fastq ;
    """
}
