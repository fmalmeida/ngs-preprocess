process lighter {
   publishDir "${params.outdir}/illumina/after_trimming", mode: 'copy',
   saveAs: {filename ->
   // This line saves the files with specific sufixes in specific folders
            if (filename.indexOf(".cor{.fq.gz, .fq}") > 0) "reads/lighter_corrected/$filename"
            else "reads/lighter_corrected/$filename"}
   container 'fmalmeida/ngs-preprocess'
   tag "Executing Ligther (read correction) step"

   input:
    file reads
    val threads

   output:
    tuple val('lighter'), file("*_1.cor.fq.gz"), file("*_2.cor.fq.gz") optional true
    file "*.cor.fq.gz" optional true
    file 'fastqc_after_correction'

   when:
   (params.lighter_execute)

   script:
   // Check if alpha is given
   alpha_param = (params.lighter_alpha) ? "-k ${params.lighter_kmer} ${params.lighter_alpha}" : "-K ${params.lighter_kmer}"
   // Check reads library
   if (params.shortreads_type == 'paired') {
        param = "-r ${reads[1]} -r ${reads[2]}"
        quality = "*_1.cor.fq.gz *_2.cor.fq.gz"
      }
      else if (params.shortreads_type == 'single') {
        param = "-r ${reads}"
        quality = "*.cor{.fq.gz, .fq}"
      }
   """
   lighter ${param} ${alpha_param} ${params.lighter_genomeSize};
   mkdir fastqc_after_correction ;
   fastqc -t ${threads} -o fastqc_after_correction -q ${quality}
   """
}
