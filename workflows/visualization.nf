#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process merge_virus_scores {
    publishDir "$params.results/visualization", mode: 'copy', overwrite: true
    
    input:
        path score_files
        path sample_table 
    
    output:
        path "merged_scores.csv"
    
    script:
    """
    python3 - <<EOF
import pandas as pd
import glob

sample_df = pd.read_csv("${sample_table}")
sample_map = dict(enumerate(sample_df['sample_ID']))

dfs = []
file_list = "${score_files}".split()

for file in file_list:
    try:
        df = pd.read_csv(file)

        sample_num = int(file.replace('.scores.csv', ''))
        sample_name = sample_map[sample_num]
        
        if 'Species' not in df.columns:
            print(f"Warning: Species column not found in {file}")
            continue
            
        df_clean = df.set_index('Species')
        df_clean.columns = [sample_name]  # Use actual sample ID
        dfs.append(df_clean)
        
    except Exception as e:
        print(f"Error processing {file}: {str(e)}")
        continue

if dfs:
    merged_df = pd.concat(dfs, axis=1)
    merged_df = merged_df.reindex(columns=sample_df['sample_ID'])
    merged_df = merged_df.fillna(0)
    merged_df.to_csv('merged_scores.csv')
else:
    exit(1)
EOF
    """
}

process generate_virus_score_heatmap {
    publishDir "$params.results/visualization", mode: 'copy', overwrite: true
    
    input:
        path merged_scores
    
    output:
        path "virus_score_heatmap.html"   
    
    script:
    """
    python3 - <<EOF
import pandas as pd
import plotly.express as px

df = pd.read_csv('${merged_scores}')

df = df.set_index('Species')

df_all = df

row_height = 20 
plot_height = max(600, len(df_all) * row_height)

fig = px.imshow(df_all,
                labels=dict(x="Sample ID", y="Species", color="Score"),
                x=df_all.columns,
                y=df_all.index,
                color_continuous_scale='Viridis',
                aspect='auto',
                zmin=0,
                zmax=100
                )

fig.update_layout(
    title='Virus Scores Heatmap: All Species',
    height=plot_height,
    width=1800,
    yaxis=dict(tickfont=dict(size=10)),
    xaxis=dict(tickangle=35, tickfont=dict(size=8))
)

fig.write_html('virus_score_heatmap.html')
EOF
    """
}

process generate_zscore_heatmap {
    publishDir "$params.results/visualization", mode: 'copy', overwrite: true
    
    input:
        path zscore_file
        path sample_table
        path peptide_table
    
    output:
        path "zscore_heatmap.html"
    
    script:
    """
    python3 - <<EOF
import pandas as pd
import plotly.graph_objects as go
import numpy as np

zscore_files = "${zscore_file}".split()
zscore_file = [f for f in zscore_files if 'zscore' in f][0]
cpm = pd.read_csv(zscore_file, index_col=0)
sample_table = pd.read_csv("${sample_table}")
peptide_table = pd.read_csv("${peptide_table}")

sample_ids = sample_table["sample_ID"].values
cpm.columns = sample_ids  # assumes order matches
cpm = cpm.apply(pd.to_numeric, errors="coerce")

peptide_table["index"] = peptide_table.index  # preserve original index
merged = peptide_table.merge(cpm, left_index=True, right_index=True)

viruses = merged["Organism"].dropna().unique()

data_traces = []
buttons = []

for virus in viruses:
    df_virus = merged[merged["Organism"] == virus]

    if df_virus.empty:
        continue

    zscores = df_virus[sample_ids].apply(pd.to_numeric, errors='coerce')
    oligo_labels = df_virus["oligo"] if "oligo" in df_virus else df_virus["index"]
    oligo_indices = [f"{i}" for i in range(len(oligo_labels))]

    heatmap = go.Heatmap(
        z=zscores.values,
        x=sample_ids,
        y=oligo_indices,
        coloraxis="coloraxis",
        visible=False,
        name=virus,
        text=oligo_labels,
        hovertemplate="Sample: %{x}<br>Peptide: %{text}<br>Z-score: %{z}<extra></extra>"
    )
    data_traces.append(heatmap)

    buttons.append(dict(
        label=virus,
        method="update",
        args=[{"visible": [i == len(data_traces) - 1 for i in range(len(data_traces))]},
              {"title": f"Z-score Heatmap: {virus}"}]
    ))

if data_traces:
    data_traces[0].visible = True
    buttons[0]["args"][0]["visible"][0] = True

    layout = go.Layout(
        title=f"z-scores Heatmap: {viruses[0]}",
        updatemenus=[dict(
            active=0,
            buttons=buttons,
            direction="down",
            showactive=True,
            x=0.5,
            xanchor="center",
            y=1.1,
            yanchor="bottom"
        )],
        coloraxis=dict(
            colorscale="Viridis",
            colorbar=dict(title="Z-score"),
        ),
        yaxis=dict(title="oligo", tickfont=dict(size=12)),
        xaxis=dict(title="Sample", tickangle=45),
        height=1850,
        annotations=[
            dict(
                text="B = Baseline E = Endline neg1 & neg2 = Negative Controls",
                xref="paper",
                yref="paper",
                x=0.5,
                y=-0.18,
                showarrow=False,
                font=dict(size=12),
                xanchor="center"
            )
        ]
    )

    fig = go.Figure(data=data_traces, layout=layout)
    fig.write_html('zscore_heatmap.html')
else:
    exit(1)
EOF
    """
}

workflow VISUALIZE {
    take:
        virus_scores
        zscore_file
        sample_table
        peptide_table
    
    main:
        merge_virus_scores(virus_scores, sample_table)  
        generate_virus_score_heatmap(merge_virus_scores.out)
        generate_zscore_heatmap(zscore_file, sample_table, peptide_table)
    
    emit:
        virus_score_heatmap = generate_virus_score_heatmap.out
        zscore_heatmap = generate_zscore_heatmap.out
}