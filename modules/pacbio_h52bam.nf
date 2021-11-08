process h52bam {
  publishDir "${params.outdir}/longreads/${id}", mode: 'copy'
  tag "Converting pacbio legacy h5 files to bam"

  input:
  file h5bas

  output:
  file("*.subreads.bam") // Get all bam files produced

  script:
  id = file("*.subreads.bam").getBaseName() - ".subreads.bam"
  """  
  # Produce bam
  bax2bam ${h5bas}/*.bas.h5 --subread --allowUnrecognizedChemistryTriple --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,SubstitutionQV,PulseWidth,SubstitutionTag;
  """
}
