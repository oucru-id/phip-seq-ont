#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process zscore_fdr_analysis {
    publishDir "$params.results/wide_data/", mode: 'copy', overwrite: true
    
    input:
        path(zscore_file)

    output:
        path "*.csv.gz"
    
    when:
        params.run_zscore_fdr_analysis
    
    script:
    """
    ZSCORE_FILE="\$(basename ${zscore_file} .gz)"
    
    if [[ "${zscore_file}" == *.gz ]]; then
        gunzip -c "${zscore_file}" > "\$ZSCORE_FILE"
    else
        cp "${zscore_file}" "\$ZSCORE_FILE"
    fi

    python3 - <<EOF
import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

z_scores = pd.read_csv("\$ZSCORE_FILE", index_col=0)
p_values = z_scores.copy()
p_values[:] = 2 * (1 - norm.cdf(np.abs(z_scores)))
adjusted_pvals = pd.DataFrame(index=z_scores.index, columns=z_scores.columns)
for sample in p_values.columns:
    _, adj_p, _, _ = multipletests(p_values[sample].values, method='fdr_bh')
    adjusted_pvals[sample] = adj_p
significant_hits = adjusted_pvals < 0.05
oligo_hit_counts = significant_hits.sum(axis=1)
sample_hit_counts = significant_hits.sum(axis=0)

adjusted_pvals.to_csv("adjusted_pvalues.csv")
log_adj_pvals = -np.log10(adjusted_pvals.astype(float))
log_adj_pvals.replace([np.inf, -np.inf], np.nan, inplace=True)
log_adj_pvals.fillna(0, inplace=True)
log_adj_pvals.to_csv("log_adjusted_pvalues.csv")
significant_hits.to_csv("significant_hits.csv")
oligo_hit_counts.to_csv("oligo_hit_counts.csv", header=["hit_count"])
sample_hit_counts.to_csv("sample_hit_counts.csv", header=["hit_count"])
EOF

    gzip -f *.csv
    """
}

workflow FDR {
    take: 
        dump_wide_csv

    main:
        zscore_fdr_analysis(
            dump_wide_csv
                .flatten()
                .filter { it.toString().contains('zscore') }
        )

    emit:
        zscore_fdr_analysis.out
}