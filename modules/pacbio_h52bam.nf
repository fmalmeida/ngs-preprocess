process h52bam {
  publishDir "${params.outdir}/longreads", mode: 'copy'
  tag "Converting pacbio legacy h5 files to bam"

  input:
  file h5bas

  output:
  file("*/*.subreads.bam") // Get all bam files produced

  when:
  !(h5bas =~ /input.*/)

  script:
  """
  # get name
  bas_file=\$(basename ${h5bas}/*.bas.h5 .bas.h5)

  # create dir with id
  mkdir "\$bas_file"

  # Produce bam
  bax2bam ${h5bas}/*.bas.h5 --subread --allowUnrecognizedChemistryTriple --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,SubstitutionQV,PulseWidth,SubstitutionTag;

  # save in directory
  mv *.subreads.bam \$bas_file
  """
}
