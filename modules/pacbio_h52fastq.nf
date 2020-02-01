process pacbio_h52fastq {
  publishDir "${params.outdir}/pacbio", mode: 'copy'
  container 'fmalmeida/ngs-preprocess'
  tag "Extracting FASTQ from pacbio legacy h5 files"

  input:
    file h5bas
    file h5bax
  output:
    file "${id}.fastq"

  script:
  id = (h5bas.getBaseName() - ".bas")
  """
  bash5tools.py --outFilePrefix ${id} --readType subreads \
  --outType fastq --minLength 200 ${h5bas} ;
  """
}
