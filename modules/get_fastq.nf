process GET_FASTQ {
  publishDir "${params.output}/SRA_FETCH/FASTQ", mode: 'copy'
  tag "${sra_ids.replaceAll(~/\s/,'')}"
  label 'process_low'
  
  input:
  val(sra_ids)

  output:
  path "*_data", emit: sra_fastq

  when: sra_ids

  script:
  def sra_ids = "${sra_ids.replaceAll(~/\s/,'')}"
  """
  fasterq-dump \\
    --split-files \\
    --threads $task.cpus \\
    --outdir ./${sra_ids}_data \\
    --progress \\
    $sra_ids &> fasterq-dump.err || true
  
  # is pacbio?
  if [ \$(grep -ic "is PACBIO, please use fastq-dump instead" fasterq-dump.err) -eq 1 ]
  then
      fastq-dump \\
        --gzip \\
        --outdir ./${sra_ids}_data \\
        $sra_ids
  else
    echo "fasterq-dump error was:"
    cat fasterq-dump.err
  fi

  # make sure they have right extension
  for i in \$( find ./${sra_ids}_data -name "*.fq" ) ; do
    mv \$i \${i%%.fq}.fastq ;
  done
  for i in \$( find ./${sra_ids}_data -name "*.fq.gz" ) ; do
    mv \$i \${i%%.fq.gz}.fastq.gz ;
  done

  # make sure data is compressed
  for i in \$( find ./${sra_ids}_data -name "*.fastq" ) ; do
    gzip \$i ;
  done
  """
}
