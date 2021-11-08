process flash {
  publishDir "${params.outdir}/shortreads/${id}/after_trimming", mode: 'copy',
       saveAs: {filename ->
  // This line saves the files with specific sufixes in specific folders
         if (filename.indexOf(".fastq") > 0) "flash_merged/$filename"
         else "flash_merged/$filename" }
  tag "Executing FLASH read merger"

  input:
  tuple val(id), file(read1), file(read2)
  
  output:
  file "flash_merged*"

  when:
  !(read1 =~ /input.*/) && !(read2 =~ /input.*/)

  script:
  """
  # run FLASH
  flash -q -o flash_merged -z -t ${params.threads} ${read1} ${read2} &> flash.log;
  """
}
