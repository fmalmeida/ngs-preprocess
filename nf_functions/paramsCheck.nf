
def paramsCheck() {
    // check lighter genome size
    if (params.lighter && !params.lighter_genome_size) {
      println """
      ERROR!
      A major error has occurred!
        ==> When using lighter, the expected genome size must be provided. Please provide an
        expected genome size with --lighter_genome_size.
      Cheers.
      """.stripIndent()
      exit 1
    }
}
