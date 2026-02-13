#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process validate_sample_table {
    input: path samples
    output: path "validated_sample_table.csv"
    script:
    """
    validate-sample-table.py \
        -s $samples \
        -o validated_sample_table.csv \
        --run_zscore_fit_predict ${params.run_zscore_fit_predict}
    """  
}

process validate_peptide_table{
    input: path peptides
    output: path "validated_peptide_table.csv"
    script:
    """
    validate-peptide-table.py \
        -p $peptides \
        -o validated_peptide_table.csv
    """
}

process generate_fasta_reference {
    input: path peptide_table
    output: path "peptides.fasta"
    script:
    """
    generate-fasta.py \
        -pt $peptide_table \
        -o peptides.fasta
    """
}

process generate_index {
    input:
    path oligo_fasta
    output:
    tuple val("peptide_ref"), path("peptide_index")
    shell:
    """
    mkdir -p peptide_index
    bowtie2-build \
        --threads ${task.cpus} \
        ${oligo_fasta} \
        peptide_index/peptide
    """
}

process nanopore_alignment {
    label 'alignment_tool'
    input:
    tuple val(sample_id), path(index), path(respective_replicate_path)
    output:
    tuple val(sample_id), path("${sample_id}.sam")
    shell:
    template "nanopore_alignment.sh"

}

process sam_to_stats {
    input:
    tuple val(sample_id), path(sam_file)
    output:
    path "${sample_id}.stats"
    shell:
    template "sam_to_stats.sh"
}


process sam_to_counts {
    input: tuple val(sample_id), path(sam_file)
    output: path "${sample_id}.counts"
    shell:
    template "sam_to_counts.sh"
}


process collect_phip_data {
    input:
    path counts_files
    path stats_files 
    path sample_table
    path peptide_table

    output:
    path "data.phip"

    script:
    """
    merge-counts-stats.py \
        -st ${sample_table} \
        -pt ${peptide_table} \
        -cfp "*.counts" \
        -sfp "*.stats" \
        -o data.phip \
    """
}

process replicate_counts {
    input: path ds
    output: path "replicated_counts.phip"
    script: 
    """
    replicate-counts.py \
        -ds ${ds} \
        -o replicated_counts.phip
    """
}

workflow ALIGN {

    main:
        sample_ch = Channel.fromPath(params.sample_table)
        peptide_ch = Channel.fromPath(params.peptide_table)

        validate_sample_table(sample_ch)
        validate_peptide_table(peptide_ch) \
            | generate_fasta_reference | generate_index

        validate_sample_table.out
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    "peptide_ref",
                    row.sample_id,
                    file(
                        row.fastq_filepath.startsWith('/') ? row.fastq_filepath : "$params.reads_prefix/$row.fastq_filepath",
                        checkIfExists: true
                    )
                )
            }
            .set { samples_ch }

        nanopore_alignment(
            generate_index.out
                .cross(samples_ch)
                .map { ref, sample ->
                    tuple(
                        sample[1],                    
                        ref[1],                      
                        file(sample[2])              
                    )
                }
        ) | (sam_to_counts & sam_to_stats)

        ds = collect_phip_data(
            sam_to_counts.out.collect(),      
            sam_to_stats.out.collect(),       
            validate_sample_table.out,
            validate_peptide_table.out
        )


        final_output = ds
        if ( params.replicate_sequence_counts )
            final_output = replicate_counts(ds)

    emit:
        final_output
}
