process pacbio_h52fastq {
  publishDir "${params.outdir}/longreads/pacbio", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Extracting FASTQ from pacbio legacy h5 files"

  input:
    file h5bas

  output:
    file "${id}.fastq" // Save fastq
    file("*.subreads.bam") // Get all bam files produced

  script:
  id = (h5bas.getBaseName() - ".bas")
  """
  bash5tools.py --outFilePrefix ${id} --readType subreads \
  --outType fastq --minLength 200 ${h5bas} ;

  # Also produce bam
  source activate pbtools ;
  bax2bam ${h5bas} --subread --allowUnrecognizedChemistryTriple \
  --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,SubstitutionQV,PulseWidth,SubstitutionTag;
  """
}
