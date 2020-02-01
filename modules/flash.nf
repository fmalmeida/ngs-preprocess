process flash {
  publishDir "${params.outdir}/shortreads/after_trimming", mode: 'copy',
       saveAs: {filename ->
  // This line saves the files with specific sufixes in specific folders
         if (filename.indexOf(".fastq") > 0) "reads/flash_merged/$filename"
         else "reads/flash_merged/$filename" }
  container 'fmalmeida/ngs-preprocess'
  tag "Executing FLASH read merger"

  input:
    file reads
    val threads
  output:
    file "flash_merge*"

  when:
  (params.flash_execute)

  script:
  """
  source activate flash ;
  flash -q -o flash_merge -z -t ${threads} ${reads[1]} ${reads[2]} &> flash.log;
  """
}
