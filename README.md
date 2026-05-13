# PhIP-Seq Analysis pipeline from ONT reads

A Nextflow pipeline for Phage Immunoprecipitation Sequencing (PhIP-Seq) analysis from Oxford Nanopore Technologies (ONT) reads, producing enriched epitopes signal, pathogen exposure scores, IEDB epitope annotations, HL7 FHIR R4 reports, and interactive 3D protein structure visualizations. [Full documentation](https://phipseq-pipeline-docs.readthedocs.io/)

## Installation

**From the repo**

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Clone the repository
git clone https://github.com/oucru-id/phip-seq-ont.git
cd phip-seq-ont

# Verify Nextflow installation
nextflow -v
```

## Directory Structure

```
phip-seq-ont
├── main.nf                          # Main workflow entry point
├── nextflow.config                  # Configuration and parameters
├── workflows/
│   ├── alignment.nf                 # Read alignment and count collection
│   ├── statistics.nf                # Normalization and Z-score calculation
│   ├── output.nf                    # dataset export
│   ├── fhir_report.nf               # HL7 FHIR R4 bundle generation
│   ├── fdr.nf                       # Benjamini-Hochberg correction
│   ├── virusscore.nf                # Per-species exposure scoring
│   ├── iedb_annotation.nf           # IEDB Epitope annotation and novelty classification
│   ├── visualization.nf             # Interactive heatmaps
│   ├── streamlit.nf                 # 3D protein visualization dashboard
│   ├── neutralization_score.nf      # Neutralization scoring (in development)
│   ├── aggregate.nf                 # result aggregation
│   └── edgeR_BEER.nf                # Optional edgeR/BEER analysis
├── bin/
│   ├── calc_scores_nofilter.py      # Virus score calculation
│   ├── fit-predict-zscore.py        # Z-score calculation
│   ├── generate-fasta.py            # Peptide FASTA generation
│   ├── merge-counts-stats.py        # Count/stats merging into phippery dataset
│   ├── phipseq_to_fhir.py           # PhIP-Seq to FHIR converter
│   ├── replicate-counts.py          # Replicate sequence count aggregation
│   ├── run_BEER.Rscript             # BEER Bayesian enrichment (optional)
│   ├── run_edgeR.Rscript            # edgeR differential enrichment (optional)
│   ├── validate-peptide-table.py    # Peptide table validation
│   └── validate-sample-table.py    # Sample table validation
├── templates/                       # Shell script templates (alignment, SAM processing)
└── data/                            # dataset (sample table, peptide table, FASTQs)
```

## Input Data

### Sample Table

A CSV file with one row per sample

| Column | Description |
|---|---|
| `sample_ID` | Sample identifier |
| `fastq_filepath` | Path to the sample's FASTQ file |
| `control_status` | `empirical`, `beads_only` (≥2 required), or `library` |

### Peptide Table

A CSV file describing the phage display library.

| Column | Description |
|---|---|---|
| `oligo` | Oligonucleotide sequence |
| `Organism` | Virus/organism name |
| `Species` | Species-level grouping |
| `peptide` | Amino acid sequence |
| `pdb_id` | PDB identifier for 3D visualization and neutralizing scoring |

### Sequencing Data

- Format: FASTQ (`.fastq.gz` or `.fastq`)
- Platform: Oxford Nanopore Technologies (ONT)
- Set `params.fastq_stream_func = 'zcat'` for gzipped files or `'cat'` for uncompressed

### Other Input Files

| File | Parameter | Description |
|---|---|---|
| IEDB database | `params.iedb_database` | TSV from [iedb.org](https://www.iedb.org/) for epitope annotation |
| PDB structures | `params.pdb_dir` | Directory of `.pdb`/`.cif` files for 3D visualization and neutralizing scoring |
| Neutralization DBs | `params.neutralization_sars_db` etc. | TSV files for neutralization scoring |

## Usage

### Run with Example Data

```bash
nextflow run main.nf
```

### Run with Custom Data

```bash
nextflow run main.nf \
  --sample_table /path/to/sample_table.csv \
  --peptide_table /path/to/peptide_table.csv \
  --results /path/to/results/
```

### Run with Docker

```bash
nextflow run main.nf -profile docker
```

### Deploy Streamlit Dashboard

After the pipeline completes:

```bash
cd results/streamlit_app/
chmod +x deploy_streamlit.sh
./deploy_streamlit.sh
```

## Support

[GitHub Issues](https://github.com/oucru-id/phip-seq-ont/issues)
