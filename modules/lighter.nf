process lighter {
   publishDir "${params.outdir}/shortreads/after_trimming", mode: 'copy',
   saveAs: {filename ->
   // This line saves the files with specific sufixes in specific folders
            if (filename.indexOf(".cor{.fq.gz, .fq}") > 0) "reads/lighter_corrected/$filename"
            else "reads/lighter_corrected/$filename"}
   tag "Executing Ligther (read correction) step"

   input:
    file reads

   output:
    tuple val('lighter'), file("*_1.cor.fq.gz"), file("*_2.cor.fq.gz") optional true
    file "*.cor.fq.gz" optional true
    file 'fastqc_after_correction'

   script:
   // Check if alpha is given
   alpha_param = (params.lighter_alpha) ? "-k ${params.lighter_kmer} ${params.lighter_genome_size} ${params.lighter_alpha}" : "-K ${params.lighter_kmer} ${params.lighter_genome_size}"
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
   # run lighter
   lighter ${param} ${alpha_param};
   mkdir fastqc_after_correction ;
   fastqc -t ${params.threads} -o fastqc_after_correction -q ${quality}
   """
}
