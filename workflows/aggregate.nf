#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


process split_samples {

    input:
        path "*"

    output: path "sample_list"
    when: params.summarize_by_organism
    shell:
    template "split_samples.py"
}

process aggregate_organisms {
    tag "${sample_id}"
    cpus 1
    memory "4.GB"
    input:
        tuple path("*"), val(sample_id)
        path public_epitopes_csv
    output: path "*.csv.gz"
    when: params.summarize_by_organism
    shell:
    template "aggregate_organisms.py"
}

process join_organisms {
    publishDir "$params.results/aggregated_data/", mode: 'copy', overwrite: true
    input: path "input/"
    output: path "*.csv.gz"
    when: params.summarize_by_organism
    shell:
    template 'join_organisms.py'
}

workflow AGG {
    take:
        dump_binary
        dump_wide_csv
        dump_tall_csv
    main:

    split_samples(dump_wide_csv)

    aggregate_organisms(
        dump_wide_csv
            .toSortedList()
            .combine(
                split_samples
                    .out
                    .splitText(){it.replace("\n", "")}
            ),
        file("${params.public_epitopes_csv}")
    )

    join_organisms(
        aggregate_organisms
            .out
            .flatten()
            .toSortedList()
    )
}


