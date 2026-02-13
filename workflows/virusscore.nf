#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process split_significant_hits {
    publishDir "$params.results/virusscore/split_hits", mode: 'copy', overwrite: true

    input:
        path significant_hits

    output:
        path "hits/*.csv.gz"

    when:
        params.run_virus_score_analysis

    script:
    """
    mkdir -p hits

    python3 - <<EOF
import pandas as pd
import gzip
import os

if "${significant_hits}".endswith('.gz'):
    df = pd.read_csv(gzip.open("${significant_hits}"), index_col=0)
else:
    df = pd.read_csv("${significant_hits}", index_col=0)

df.index.name = "id"

for col in df.columns:
    single_col = df[[col]].copy()
    single_col = single_col.reset_index()  
    output_file = f"hits/{col}.csv"
    single_col.to_csv(output_file, index=False)
EOF

    for f in hits/*.csv; do
        [ -f "\$f" ] && gzip -f "\$f"
    done
    """
}

process calculate_virus_scores {
    publishDir "$params.results/virusscore/scores", mode: 'copy', overwrite: true
    
    input:
        tuple path(hits_file), path(oligo_metadata), val(level), val(epitope_len)
    
    output:
        path "*.scores.csv"

    when:
        params.run_virus_score_analysis
    
    script:
    """
    SAMPLE=\$(basename ${hits_file} .csv.gz)
    
    python3 $baseDir/bin/calc_scores_nofilter.py \\
        "${hits_file}" \\
        "${oligo_metadata}" \\
        "${level}" \\
        ${epitope_len} > "\${SAMPLE}.scores.csv"
    """
}

workflow VIRUSSCORE {
    take:
        dump_wide_csv
        oligo_metadata

    main:
        significant_hits = dump_wide_csv
            .flatten()
            .filter { it.toString().contains('significant_hits') }

        split_significant_hits(significant_hits)
            .flatten()
            .set { split_files }

        split_files
            .combine(Channel.fromPath(oligo_metadata))
            .combine(Channel.value(params.virus_score_level))
            .combine(Channel.value(params.virus_score_epitope_len))
            .set { virus_score_input }

        calculate_virus_scores(
            virus_score_input
        )

    emit:
        calculate_virus_scores.out.collect()
}