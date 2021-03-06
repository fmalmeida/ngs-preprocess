process pacbio_h52bam {
  publishDir "${params.outdir}/longreads/pacbio", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Converting pacbio legacy h5 files to bam"

  input:
    file h5bas
    val h5bas_dir

  output:
    file("*.subreads.bam") // Get all bam files produced

  script:
  """
  # Produce bam
  source activate pbtools-h5 ;
  bax2bam ${h5bas_dir}/*.bas.h5 --subread --allowUnrecognizedChemistryTriple \
  --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,SubstitutionQV,PulseWidth,SubstitutionTag;
  """
}
