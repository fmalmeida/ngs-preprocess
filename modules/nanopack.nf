process NANOPACK {
    publishDir "${params.output}/preprocessing_outputs/${meta.longreads_type}/QC/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    
    output:
    path "NanoPlot"
    path "NanoQC"
    path "NanoStats"

    script:
    """
    # Plotting
    NanoPlot \\
        -t ${task.cpus} \\
        --fastq ${reads} \\
        -o NanoPlot \\
        --N50 \\
        --title "${meta.id} sample" \\
        --plots hex dot kde ;

    # Checking Quality
    nanoQC \\
        -o NanoQC \\
        ${reads} ;

    # Generate Statistics Summary
    NanoStat \\
        --fastq ${reads} \\
        -t ${task.cpus} \\
        -n ${meta.id}.txt \\
        --outdir NanoStats ;
    """
}
