process H52BAM {
  publishDir "${params.output}/preprocessing_outputs/pacbio/h52bam", mode: 'copy'
  tag "${id}"

  input:
  file h5bas

  output:
  path("*/*.subreads.bam") // Get all bam files produced
  path("*")

  when:
  !(h5bas =~ /input.*/)

  script:
  id = (h5bas.getBaseName())
  """
  # get name
  bas_file=\$(basename ${h5bas}/*.bas.h5 .bas.h5)

  # create dir with id
  mkdir "\$bas_file"

  # Produce bam
  bax2bam \\
      ${h5bas}/*.bas.h5 \\
      --subread \\
      --allowUnrecognizedChemistryTriple \\
      --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,SubstitutionQV,PulseWidth,SubstitutionTag;

  # save in directory
  mv *.subreads.bam \$bas_file
  """
}
