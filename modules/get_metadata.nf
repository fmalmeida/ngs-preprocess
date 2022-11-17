process GET_METADATA {
  publishDir "${params.output}/SRA_FETCH", mode: 'copy'
  tag "$id"
  label 'process_low'
  
  input:
  path(sra_run)
  val(headers)

  output:
  path "*.csv", emit: sra_metadata

  script:
  id = "${sra_run.baseName}" - "_data"
  file_string = "${sra_run.toRealPath().toUriString()}"
  """
  esearch -db sra -query \'${id}\' | \\
    efetch -format runinfo | \\
    csvtk cut -f ${headers} \\
    > tmp.txt
  
  paste -d',' <(sed '1q;d' tmp.txt) <(echo fastq_dir) >  ${id}_sra_runInfo.csv
  paste -d',' <(sed '2q;d' tmp.txt) <(echo ${file_string}) >> ${id}_sra_runInfo.csv
  """
}
