process flash {
  publishDir "${params.outdir}/shortreads/${id}/after_trimming", mode: 'copy',
       saveAs: {filename ->
  // This line saves the files with specific sufixes in specific folders
         if (filename.indexOf(".fastq") > 0) "flash_merged/$filename"
         else "flash_merged/$filename" }
  tag "Executing FLASH read merger"

  input:
    file reads
  output:
    file "flash_merged*"

  script:
  id = reads[0]
  """
  # run FLASH
  flash -q -o flash_merged -z -t ${params.threads} ${reads[1]} ${reads[2]} &> flash.log;
  """
}
