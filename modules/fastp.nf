process FASTP {
    publishDir "${params.output}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".fq.gz")) "final_output/short_reads/$filename"
        else if (filename.endsWith(".json")) "preprocessing_outputs/short_reads/$filename"
        else if (filename.endsWith(".html")) "preprocessing_outputs/short_reads/$filename"
        else "$filename"
    }
    tag "${meta.id}"
    label 'process_low'
    
    input:
    tuple val(meta), path(reads)

    output:
    path "*"

    script:
    if (meta.shortreads_type == 'paired') {
        if (params.fastp_merge_pairs) {
            reads_param = "--in1 ${reads[0]} --in2 ${reads[1]} --out1 ${meta.id}_R1.unmerged.fq.gz --out2 ${meta.id}_R2.unmerged.fq.gz --detect_adapter_for_pe --merge --merged_out ${meta.id}.merged.fq.gz"
        } else {
            reads_param = "--in1 ${reads[0]} --in2 ${reads[1]} --out1 ${meta.id}_R1.preprocessed.fq.gz --out2 ${meta.id}_R2.preprocessed.fq.gz --detect_adapter_for_pe"
        }
    } else {
        reads_param = "-i ${reads} -o ${meta.id}.preprocessed.fq.gz"
    }
    correction_param = (params.fastp_correct_pairs)         ? "--correction"                          : ""
    additional_param = (params.fastp_additional_parameters) ? "${params.fastp_additional_parameters}" : ""
    """
    # run fastp
    fastp \\
        $additional_param \\
        --thread ${task.cpus} \\
        --average_qual ${params.fastp_average_quality} \\
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html \\
        ${reads_param} \\
        ${correction_param}
    """
}
