process porechop {
  publishDir "${params.outdir}/longreads/porechop_out", mode: 'copy'
  tag "Trimming longreads with Porechop"

  input:
    file reads
    val threads
    val barcode
  output:
    file "${id}_trimmed.fastq" optional true
    file "porechop_barcodes/*.fastq" optional true

  script:
    id = (reads.getBaseName() - "fastq.gz" - ".fastq")
    if (barcode)
    """
    porechop -i ${reads} -b porechop_barcodes --barcode_threshold 85
    """
    else
    """
    porechop -i ${reads} -t ${threads} --format fastq -o ${id}_trimmed.fastq ;
    """
}
