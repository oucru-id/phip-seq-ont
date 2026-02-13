#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CREATE_FHIR {
    publishDir "${params.results}/fhir", mode: 'copy'
    
    input:
    path(zscore_data)
    path(sample_table)

    output:
    path "*.fhir.json", emit: fhir_output
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/phipseq_to_fhir.py \\
        --zscore_input ${zscore_data} \\
        --sample_table ${sample_table} \\
        --output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "fhir_converter":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}

workflow FHIR {
    take:
        dump_wide_csv

    main:
        sample_table_ch = Channel.fromPath(params.sample_table, checkIfExists: true)
        
        CREATE_FHIR(
            dump_wide_csv
                .flatten()
                .filter { it.toString().contains('zscore') },
            sample_table_ch
        )

    emit:
        fhir_output = CREATE_FHIR.out.fhir_output
        versions = CREATE_FHIR.out.versions
}