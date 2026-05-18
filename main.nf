nextflow.enable.dsl = 2

if (params.sample_table != "$baseDir/data/VIR3-example/sample_table_with_beads_no_lib.csv")
    params.reads_prefix = "$launchDir"
else
    params.reads_prefix = "$baseDir"

log.info """\
PhIP-Seq nanopore reads analysis pipeline
Developed by SPHERES OUCRU-ID 
Adapted from phip-flow (https://github.com/matsengrp/phip-flow)
Version: ${params.version}
"""

include { ALIGN } from './workflows/alignment.nf'
include { STATS } from './workflows/statistics.nf'
include { DSOUT } from './workflows/output.nf'
include { FHIR } from './workflows/fhir_report.nf'
include { FDR } from './workflows/fdr.nf'
include { VIRUSSCORE } from './workflows/virusscore.nf'
include { IEDB } from './workflows/iedb_annotation.nf'
include { VISUALIZE } from './workflows/visualization.nf'
include { STREAMLIT } from './workflows/streamlit.nf'
include { AGG } from './workflows/aggregate.nf'
include { NEUTRALIZATION_PREDICTION } from './workflows/neutralization_score.nf'

workflow {
    ALIGN()
    STATS(ALIGN.out)
    DSOUT(STATS.out)
    FDR(DSOUT.out[1]) 
    VIRUSSCORE(
        FDR.out,
        params.oligo_metadata
    )
    FHIR(
        DSOUT.out[1],
        file(params.peptide_table),
        VIRUSSCORE.out
    )
    
    IEDB(
        DSOUT.out[1].flatten().filter { it.toString().contains('zscore') },
        file(params.peptide_table)
    )
    
    VISUALIZE(
        VIRUSSCORE.out,  
        DSOUT.out[1].flatten().filter { it.toString().contains('zscore') },
        file(params.sample_table),
        file(params.peptide_table)
    )
    
    STREAMLIT(
        DSOUT.out[1].flatten().filter { it.toString().contains('zscore') },
        file(params.peptide_table),
        file(params.pdb_dir)
    )
    
    NEUTRALIZATION_PREDICTION(
        IEDB.out.significant,
        IEDB.out.annotated,
        file(params.peptide_table),
        DSOUT.out[1].flatten().filter { it.toString().contains('zscore') },
        file(params.pdb_dir),
        file(params.neutralization_db)
    )

    STREAMLIT.out.streamlit_app.view { 
        if (params.deploy_streamlit) {
           println """            
            Streamlit Dashboard: ${params.results}/streamlit_app/
                        
            cd ${params.results}/streamlit_app/
            ./deploy_streamlit.sh
            
            """
        }
    }
}
