# NGS-preprocess pipeline tests

## Help message

* Command: `nextflow run fmalmeida/NGS-preprocess --help`
	* Working

## Error messages

Started to set error messages

## Usages

* Pacbio input (.subreads.bam)
	* One at a time: Working
	* Multiple at once: Working.
	* Barcoded: Working
	* OBS: files .subreads.bam and .subreads.bam.pbi must ALWAYS be in the same directory
* Pacbio legacy input (*.bas.h5)
	* Working
	* Barcoded: Not parsed
* Pacbio input (.fastq)
	* Working
* ONT input (.fastq)
	* Working
* Illumina input (single-end)
	* A single input: Working
	* Multiple inputs: Working
* Illumina input (paired-end)
	* A single input pair: Working

## Report

The parameter usage for report generation is working
