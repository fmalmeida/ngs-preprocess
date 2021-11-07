process porechop {
  publishDir "${params.outdir}/longreads/porechop_out", mode: 'copy'
  tag "Trimming longreads with Porechop"

  input:
    file reads
  output:
    file "${id}_trimmed.fastq" optional true
    file "porechop_barcodes/*.fastq" optional true

  script:
    id = (reads.getBaseName() - "fastq.gz" - ".fastq")
    if (params.nanopore_is_barcoded)
    """
    porechop -i ${reads} -t ${params.threads} -b porechop_barcodes --barcode_threshold 85
    """
    else
    """
    porechop -i ${reads} -t ${params.threads} --format fastq -o ${id}_trimmed.fastq ;
    """
}
