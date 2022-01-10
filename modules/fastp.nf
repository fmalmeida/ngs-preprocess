process FASTP {
  publishDir "${params.output}", mode: 'copy', saveAs: { filename ->
    if (filename.endsWith(".fq.gz")) "final_output/illumina/$filename"
    else if (filename.endsWith(".json")) "preprocessing_outputs/illumina/QC/$filename"
    else if (filename.endsWith(".html")) "preprocessing_outputs/illumina/QC/$filename"
    else "$filename"
  }
  tag "${id}"
  
  input:
  tuple val(id), file(read1), file(read2)
  file(sreads)

  output:
  path "*"

  when:
  (!(read1 =~ /input.*/) && !(read2 =~ /input.*/)) || !(sreads =~ /input.*/)

  script:
  if (params.shortreads_type == 'paired') {
      if (params.fastp_merge_pairs) {
          reads_param = "--in1 ${read1} --in2 ${read2} --out1 ${id}_R1.unmerged.fq.gz --out2 ${id}_R2.unmerged.fq.gz --detect_adapter_for_pe --merge --merged_out ${id}.merged.fq.gz"
      } else {
          reads_param = "--in1 ${read1} --in2 ${read2} --out1 ${id}_R1.preprocessed.fq.gz --out2 ${id}_R2.preprocessed.fq.gz --detect_adapter_for_pe"
      }
  } else if (params.shortreads_type == 'single') {
      id = sreads.getBaseName() - ".fastq.gz" - ".fastq"
      reads_param = "-i ${sreads} -o ${id}.preprocessed.fq.gz"
  }
  correction_param = (params.fastp_correct_pairs) ? "--correction" : ""
  """
  # run fastp
  fastp \\
      --thread ${params.threads} \\
      --average_qual ${params.fastp_average_quality} \\
      --json ${id}_fastp.json \\
      --html ${id}_fastp.html \\
      ${reads_param} \\
      ${correction_param} \\
      ${params.fastp_additional_parameters}
  """
}