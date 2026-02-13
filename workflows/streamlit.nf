#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GENERATE_JSMOL_HTML {
    publishDir "${params.results}/interactive_viz", mode: 'copy'
    
    input:
    path zscore_file
    path peptide_info
    path pdb_dir
    
    output:
    path "*.html", emit: html_files
    path "js/", emit: js_files, optional: true
    path "pdb_files/", emit: pdb_files
    path "*.log", emit: logs
    path "epitope_data.csv", emit: epitope_data
    
    script:
    """
#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import json
import shutil
import logging
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('jsmol_generation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

try:
    os.makedirs("js", exist_ok=True)
    os.makedirs("pdb_files", exist_ok=True)
    
    logging.info("Reading zscore data from '${zscore_file}'")
    zscore_data = pd.read_csv("${zscore_file}")
    
    logging.info("Reading peptide info from '${peptide_info}'")
    peptide_info = pd.read_csv("${peptide_info}")
    
    if 'Unnamed: 0' in zscore_data.columns:
        zscore_data = zscore_data.rename(columns={'Unnamed: 0': 'id'})
    
    peptide_info = peptide_info.dropna(subset=['pdb_id'])
    peptide_info = peptide_info[peptide_info['pdb_id'].str.strip() != '']
    
    if peptide_info.empty:
        logging.error("No valid PDB IDs found")
        with open("no_structures.html", "w") as f:
            f.write("<html><body><h1>No valid PDB structures found</h1></body></html>")
        pd.DataFrame(columns=['pdb_id', 'peptide_id', 'sample', 'start', 'end', 'length', 'zscore', 'sequence']).to_csv('epitope_data.csv', index=False)
        sys.exit(0)
    
    for pdb_id in peptide_info['pdb_id'].unique():
        found = False
        for ext in ['.pdb', '.PDB', '.cif', '.CIF']:
            source_path = os.path.join("${pdb_dir}", f"{pdb_id}{ext}")
            if os.path.exists(source_path):
                target_ext = ext.lower()
                target_path = os.path.join("pdb_files", f"{pdb_id}{target_ext}")
                shutil.copy2(source_path, target_path)
                logging.info(f"Copied {source_path} to {target_path}")
                found = True
                break
        if not found:
            logging.warning(f"Structure file for {pdb_id} not found in {os.path.abspath('${pdb_dir}')}")
    
    pdb_groups = peptide_info.groupby('pdb_id')
    logging.info(f"Found {len(pdb_groups)} unique PDB structures")
    
    generated_count = 0
    all_epitope_data = []
    
    for pdb_id, group in pdb_groups:
        logging.info(f"Processing PDB: {pdb_id}")
        
        pdb_peptides = group['id'].unique()
        pdb_zscores = zscore_data[zscore_data['id'].isin(pdb_peptides)]
        zscore_columns = [col for col in zscore_data.columns if col not in ['id']]
        
        epitopes = []
        for zscore_col in zscore_columns:
            for idx, row in pdb_zscores.iterrows():
                try:
                    zscore_val = float(row[zscore_col])
                    if zscore_val > ${params.zscore_threshold}:
                        matching_peptides = group[group['id'] == row['id']]
                        if len(matching_peptides) > 0:
                            peptide_data = matching_peptides.iloc[0]
                            start_pos = int(peptide_data['Prot_Start'])
                            
                            peptide_sequence = str(peptide_data.get('peptide', peptide_data.get('Sequence', '')))
                            if peptide_sequence and peptide_sequence != 'nan':
                                peptide_length = len(peptide_sequence)
                            else:
                                peptide_length = 35
                            
                            end_pos = start_pos + peptide_length - 1
                            
                            epitope_data = {
                                'start': start_pos,
                                'end': end_pos,
                                'zscore': round(float(zscore_val), 2)
                            }
                            epitopes.append(epitope_data)
                            
                            all_epitope_data.append({
                                'pdb_id': pdb_id,
                                'peptide_id': row['id'],
                                'sample': zscore_col,
                                'start': start_pos,
                                'end': end_pos,
                                'length': peptide_length,
                                'zscore': round(float(zscore_val), 2),
                                'sequence': peptide_sequence if peptide_sequence != 'nan' else 'Unknown'
                            })
                            
                            logging.info(f"Added epitope: {epitope_data}")
                except (ValueError, TypeError) as e:
                    logging.warning(f"Error processing zscore for {row['id']}: {e}")
                    continue
        
        logging.info(f"Total epitopes found for {pdb_id}: {len(epitopes)}")
        
        epitopes_json = json.dumps(epitopes, indent=2)
        
        html_content = f\"\"\"<!DOCTYPE html>
<html>
<head>
    <title>PhIP-seq Analysis - {pdb_id}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .epitope {{ margin: 10px 0; padding: 10px; background: #f0f0f0; border-radius: 5px; }}
    </style>
</head>
<body>
    <h1>PhIP-seq Analysis - {pdb_id}</h1>
    <p><strong>Epitopes found:</strong> {len(epitopes)}</p>
    <p><strong>Note:</strong> For interactive 3D visualization, use the Streamlit app.</p>
    <div id="epitopes">
        {''.join([f'<div class=\"epitope\">Residues {ep["start"]}-{ep["end"]} (Z-score: {ep["zscore"]})</div>' for ep in epitopes])}
    </div>
    <script>
        var epitopes = {epitopes_json};
        console.log('Epitopes data:', epitopes);
    </script>
</body>
</html>\"\"\"
        
        with open(f"{pdb_id}_interactive.html", "w") as f:
            f.write(html_content)
        
        generated_count += 1
        logging.info(f"Generated visualization for {pdb_id}")
    
    if generated_count == 0:
        logging.warning("No visualization files were generated")
        with open("no_structures.html", "w") as f:
            f.write("<html><body><h1>No structures could be processed</h1></body></html>")
    else:
        logging.info(f"Successfully generated {generated_count} visualizations")

    if all_epitope_data:
        epitope_df = pd.DataFrame(all_epitope_data)
        epitope_df.to_csv('epitope_data.csv', index=False)
        logging.info(f"Exported {len(all_epitope_data)} epitope records for Streamlit")
    else:
        pd.DataFrame(columns=['pdb_id', 'peptide_id', 'sample', 'start', 'end', 'length', 'zscore', 'sequence']).to_csv('epitope_data.csv', index=False)

except Exception as e:
    logging.error(f"Error in data processing: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
    """
}

process CREATE_STREAMLIT_APP {
    publishDir "${params.results}/streamlit_app", mode: 'copy'
    
    input:
    path html_files
    path pdb_files  
    path epitope_data
    
    output:
    path "streamlit_app.py", emit: app
    path "requirements.txt", emit: requirements
    path "config.toml", emit: config
    path "data/", emit: data_dir
    path "deploy_streamlit.sh", emit: deploy_script
    path "README_streamlit.md", emit: readme
    
    script:
    """
#!/usr/bin/env python3
import os
import shutil

os.makedirs("data", exist_ok=True)

html_files_list = [f for f in os.listdir('.') if f.endswith('.html')]
for html_file in html_files_list:
    shutil.copy2(html_file, 'data/')

if os.path.exists('pdb_files'):
    shutil.copytree('pdb_files', 'data/pdb_files', dirs_exist_ok=True)

if os.path.exists('epitope_data.csv'):
    shutil.copy2('epitope_data.csv', 'data/')

app_content = r'''import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import os
import base64
import json
import requests

st.set_page_config(
    page_title="PhIP-seq 3D Visualization",
    layout="wide"
)

AA_MAP = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def get_pdb_residue_mapping(pdb_content, fmt='pdb'):
    chains = {}
    lines = pdb_content.splitlines()
    
    if fmt == 'pdb':
        for line in lines:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                res_name = line[17:20].strip()
                chain_id = line[21]
                try:
                    res_seq = int(line[22:26])
                except ValueError:
                    continue 
                
                if chain_id not in chains:
                    chains[chain_id] = {'seq': '', 'ids': []}
                
                if res_name in AA_MAP:
                    chains[chain_id]['seq'] += AA_MAP[res_name]
                    chains[chain_id]['ids'].append(res_seq)

    elif fmt == 'cif' or fmt == 'mmcif':
        loop_indices = {}
        in_loop = False
        
        for line in lines:
            line = line.strip()
            if line.startswith('loop_'):
                in_loop = True
                loop_indices = {}
                current_col_idx = 0
                continue
            
            if in_loop and line.startswith('_atom_site.'):
                key = line.split('.')[1].strip()
                loop_indices[key] = current_col_idx
                current_col_idx += 1
                continue
            
            if in_loop and (line.startswith('ATOM') or line.startswith('HETATM')):
                if 'label_atom_id' in loop_indices and 'label_comp_id' in loop_indices and 'label_asym_id' in loop_indices and 'label_seq_id' in loop_indices:
                    parts = line.split()
                    try:
                        atom_idx = loop_indices['label_atom_id']
                        if len(parts) > atom_idx and parts[atom_idx] == 'CA':
                            res_name = parts[loop_indices['label_comp_id']]
                            chain_id = parts[loop_indices['label_asym_id']]
                            res_seq = int(parts[loop_indices['label_seq_id']])
                            
                            if chain_id not in chains:
                                chains[chain_id] = {'seq': '', 'ids': []}
                            
                            if res_name in AA_MAP:
                                chains[chain_id]['seq'] += AA_MAP[res_name]
                                chains[chain_id]['ids'].append(res_seq)
                    except (ValueError, IndexError):
                        continue
    return chains

def find_sequence_locations(peptide_seq, pdb_mapping, first_match_only=False):
    matches = []
    if not peptide_seq or peptide_seq == 'Unknown' or pd.isna(peptide_seq):
        return matches
        
    peptide_seq = str(peptide_seq).strip()
    
    found_exact = False
    for chain_id, data in pdb_mapping.items():
        seq = data['seq']
        start = 0
        while True:
            idx = seq.find(peptide_seq, start)
            if idx == -1:
                break
            
            end_idx = idx + len(peptide_seq)
            residues = data['ids'][idx:end_idx]
            matches.append({'chain': chain_id, 'resi': residues, 'method': 'exact'})
            found_exact = True
            
            if first_match_only:
                return matches
            start = idx + 1
            
    if found_exact:
        return matches

    WINDOW_SIZE = 7
    if len(peptide_seq) < WINDOW_SIZE:
        WINDOW_SIZE = len(peptide_seq)
        
    fragment_matches = {} 
    
    for i in range(len(peptide_seq) - WINDOW_SIZE + 1):
        fragment = peptide_seq[i : i + WINDOW_SIZE]
        
        for chain_id, data in pdb_mapping.items():
            seq = data['seq']
            idx = seq.find(fragment)
            if idx != -1:
                if chain_id not in fragment_matches:
                    fragment_matches[chain_id] = set()
                res_ids = data['ids'][idx : idx + WINDOW_SIZE]
                fragment_matches[chain_id].update(res_ids)

    if fragment_matches:
        for chain_id, res_set in fragment_matches.items():
            matches.append({
                'chain': chain_id, 
                'resi': sorted(list(res_set)),
                'method': 'partial'
            })
            if first_match_only:
                break
                
    return matches

st.title("PhIP-seq 3D Visualization")
st.markdown("Interactive 3D PhIP-Seq enriched epitope visualization")

@st.cache_data
def load_epitope_data():
    try:
        df = pd.read_csv('data/epitope_data.csv')
        if 'sample' in df.columns:
            df['sample'] = df['sample'].astype(str)
        return df
    except FileNotFoundError:
        return pd.DataFrame()

@st.cache_data
def fetch_pdb_content(pdb_id):
    # Try local PDB
    local_path_pdb = f'data/pdb_files/{pdb_id}.pdb'
    if os.path.exists(local_path_pdb):
        with open(local_path_pdb, 'r') as f:
            return f.read(), 'pdb'
            
    local_path_cif = f'data/pdb_files/{pdb_id}.cif'
    if os.path.exists(local_path_cif):
        with open(local_path_cif, 'r') as f:
            return f.read(), 'cif' # Return 'cif' for Mol* compatibility
    
    try:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.text, 'pdb'
        else:
            return None, None
    except Exception as e:
        st.error(f"Error fetching PDB {pdb_id}: {str(e)}")
        return None, None

st.sidebar.header("Controls")

epitope_df = load_epitope_data()

if not epitope_df.empty:
    available_pdbs = sorted(epitope_df['pdb_id'].unique())
    selected_pdb = st.sidebar.selectbox("Select Structure:", available_pdbs)
    
    pdb_subset_df = epitope_df[epitope_df['pdb_id'] == selected_pdb]
    
    if 'sample' in pdb_subset_df.columns:
        available_samples = ["View Structure Only"] + sorted(pdb_subset_df['sample'].unique())
        selected_sample = st.sidebar.selectbox("Select Sample:", available_samples)
    else:
        selected_sample = "View Structure Only"
        st.sidebar.warning("No sample information found in data.")
    
    bg_color = st.sidebar.selectbox("Background:", ["white", "black"])
    bg_color_hex = "#ffffff" if bg_color == "white" else "#000000"
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("**Visualization Options**")
    
    spin_on_load = st.sidebar.checkbox("Auto-Spin", value=False)
    
    highlight_all_chains = st.sidebar.checkbox(
        "Highlight all chains (Multimer)", 
        value=True,
        help="Highlights the sequence in all protein chains."
    )
    
    if selected_sample and selected_sample != "View Structure Only":
        filtered_df = pdb_subset_df[
            (pdb_subset_df['sample'] == selected_sample)
        ]
    else:
        filtered_df = pd.DataFrame(columns=pdb_subset_df.columns)
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Total PhIP-Seq enriched epitopes", len(filtered_df))
    with col2:
        st.metric("PDB/CIF ID", selected_pdb)
    
    st.subheader(f"3D Molecular Viewer - {selected_pdb}")
    
    with st.spinner(f'Loading PDB structure for {selected_pdb}...'):
        pdb_content, pdb_format = fetch_pdb_content(selected_pdb)
    
    if pdb_content:
        view_col, details_col = st.columns([3, 1])
        
        # Parse PDB mapping
        pdb_mapping = get_pdb_residue_mapping(pdb_content, pdb_format)
        
        selected_epitope_id = None
        
        with details_col:
            st.markdown("### Epitope Selection")
            
            if not filtered_df.empty:
                epitope_options = []
                epitope_map = {} 
                
                for _, row in filtered_df.iterrows():
                    peptide_id = str(row['peptide_id'])
                    if peptide_id.endswith('.0'):
                        peptide_id = peptide_id[:-2]
                    
                    option_label = f"Peptide ID {peptide_id} (Z={row['zscore']:.2f})"
                    epitope_options.append(option_label)
                    epitope_map[option_label] = row
                
                epitope_options.insert(0, "View All Enriched")
                
                selected_option = st.selectbox("Select Epitope:", epitope_options, index=0)
                
                if selected_option != "View All Enriched":
                    selected_row = epitope_map[selected_option]
                    selected_epitope_id = selected_row['peptide_id']
                    st.markdown("---")
                    st.markdown(f"**Peptide ID:** {selected_epitope_id}")
                    st.markdown(f"**Z-score:** {selected_row['zscore']:.2f}")
                    st.code(selected_row.get('sequence', 'N/A'))
                else:
                    st.markdown("---")
                    st.info(f"Showing all {len(filtered_df)} enriched epitopes.")
            else:
                st.markdown("*No epitopes available.*")
            
            st.markdown("---")
            with st.expander("Structure Info", expanded=False):
                for chain, data in pdb_mapping.items():
                    if data['ids']:
                        st.text(f"Chain {chain}: {min(data['ids'])} - {max(data['ids'])}")

        with view_col:
            molstar_selections = []
            
            for chain_id in pdb_mapping.keys():
                molstar_selections.append({
                    'struct_asym_id': chain_id,
                    'color': {'r': 220, 'g': 220, 'b': 220}, # Light grey
                    'focus': False
                })
            
            if not filtered_df.empty:
                if selected_epitope_id is not None:
                    rows_to_highlight = filtered_df[filtered_df['peptide_id'] == selected_epitope_id]
                    color_rgb = {'r': 0, 'g': 255, 'b': 0} 
                else:
                    rows_to_highlight = filtered_df
                    color_rgb = {'r': 255, 'g': 0, 'b': 0}
                
                for _, row in rows_to_highlight.iterrows():
                    peptide_seq = row.get('sequence', '')
                    matches = find_sequence_locations(
                        peptide_seq, 
                        pdb_mapping, 
                        first_match_only=(not highlight_all_chains)
                    )
                    
                    for match in matches:
                        residues = sorted(match['resi'])
                        if not residues: continue
                        
                        ranges = []
                        start = residues[0]
                        prev = residues[0]
                        
                        for r in residues[1:]:
                            if r != prev + 1:
                                ranges.append((start, prev))
                                start = r
                            prev = r
                        ranges.append((start, prev))
                        
                        for r_start, r_end in ranges:
                            molstar_selections.append({
                                'struct_asym_id': match['chain'],
                                'start_residue_number': r_start,
                                'end_residue_number': r_end,
                                'color': color_rgb,
                                'focus': False
                            })

            b64_pdb = base64.b64encode(pdb_content.encode()).decode()
            
            html_code = (
                f'<!DOCTYPE html>'
                f'<html lang="en">'
                f'<head>'
                f'    <meta charset="utf-8" />'
                f'    <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">'
                f'    <link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/pdbe-molstar@3.2.0/build/pdbe-molstar.css">'
                f'    <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/pdbe-molstar@3.2.0/build/pdbe-molstar-plugin.js"></script>'
                f'    <style>'
                f'        * {{ margin: 0; padding: 0; box-sizing: border-box; }}'
                f'        .viewer-section {{ position: relative; width: 100%; height: 600px; }}'
                f'        #myViewer {{ position: absolute; top: 0; left: 0; right: 0; bottom: 0; }}'
                f'    </style>'
                f'</head>'
                f'<body>'
                f'    <div class="viewer-section">'
                f'        <div id="myViewer"></div>'
                f'    </div>'
                f'    <script>'
                f'        var pdbBase64 = "{b64_pdb}";'
                f'        var pdbFormat = "{pdb_format}";'
                f'        var bgColorHex = "{bg_color_hex}";'
                f'        function hexToRgb(hex) {{'
                f'            var result = /^#?([a-f0-9]{{2}})([a-f0-9]{{2}})([a-f0-9]{{2}})\$/i.exec(hex);'
                f'            return result ? {{ r: parseInt(result[1], 16), g: parseInt(result[2], 16), b: parseInt(result[3], 16) }} : {{r: 255, g: 255, b: 255}};'
                f'        }}'
                f'        function b64toBlob(b64Data, contentType="", sliceSize=512) {{'
                f'            var byteCharacters = atob(b64Data);'
                f'            var byteArrays = [];'
                f'            for (var offset = 0; offset < byteCharacters.length; offset += sliceSize) {{'
                f'                var slice = byteCharacters.slice(offset, offset + sliceSize);'
                f'                var byteNumbers = new Array(slice.length);'
                f'                for (var i = 0; i < slice.length; i++) {{ byteNumbers[i] = slice.charCodeAt(i); }}'
                f'                var byteArray = new Uint8Array(byteNumbers);'
                f'                byteArrays.push(byteArray);'
                f'            }}'
                f'            return new Blob(byteArrays, {{type: contentType}});'
                f'        }}'
                f'        var blob = b64toBlob(pdbBase64, "text/plain");'
                f'        var blobUrl = URL.createObjectURL(blob);'
                f'        var viewerInstance = new PDBeMolstarPlugin();'
                f'        var options = {{'
                f'            customData: {{ url: blobUrl, format: pdbFormat, binary: false }},'
                f'            bgColor: hexToRgb(bgColorHex),'
                f'            hideControls: true,'
                f'            hideCanvasControls: ["selection", "controlInfo"]'
                f'        }};'
                f'        var viewerContainer = document.getElementById("myViewer");'
                f'        viewerInstance.render(viewerContainer, options);'
                f'        viewerInstance.events.loadComplete.subscribe(() => {{'
                f'            var selections = {json.dumps(molstar_selections)};'
                f'            if(selections.length > 0) {{'
                f'                try {{ viewerInstance.visual.select({{ data: selections }}); }} catch(err) {{ console.error("Selection error", err); }}'
                f'            }}'
                f'            var spin = {str(spin_on_load).lower()};'
                f'            if(spin) {{ viewerInstance.visual.toggleSpin(true); }}'
                f'        }});'
                f'    </script>'
                f'</body>'
                f'</html>'
            )
            
            components.html(html_code, height=600)
            
    else:
        st.error(f"Could not load structure for {selected_pdb}")

else:
    st.warning("No epitope data found.")

st.markdown("---")
st.markdown("*Generated by PhIP-Seq pipeline (SPHERES Lab Team) | Powered by Streamlit & MolStar | Development Phase (Prototype)*")
'''

with open('streamlit_app.py', 'w') as f:
    f.write(app_content)

with open('requirements.txt', 'w') as f:
    f.write('''streamlit>=1.28.0
pandas>=1.5.0
numpy>=1.24.0
requests>=2.28.0
''')

with open('config.toml', 'w') as f:
    f.write('''[theme]
primaryColor = "#3498db"
backgroundColor = "#ffffff"

[server]
port = 8501
maxUploadSize = 200
''')

with open('deploy_streamlit.sh', 'w') as f:
    shebang = "#!/bin/bash"
    f.write(f'''{shebang}
echo "Starting PhIP-seq Streamlit Dashboard with Mol*"

if ! command -v streamlit &> /dev/null; then
    echo "Installing dependencies"
    pip install -r requirements.txt
fi

mkdir -p .streamlit
cp config.toml .streamlit/

streamlit run streamlit_app.py --server.port=8501
''')

with open('README_streamlit.md', 'w') as f:
    f.write('''# PhIP-seq Streamlit Dashboard

Interactive dashboard for PhIP-seq epitope analysis with Mol* 3D visualization.

## Features

- **Mol*** Integration: Uses PDBe Mol* Web Component
- Interactive 3D PhIP-Seq enriched epitope visualization
- Red highlighting for all enriched epitopes
- Green highlighting for selected epitope

## Quick Start

pip install -r requirements.txt
chmod +x deploy_streamlit.sh
./deploy_streamlit.sh

''')

os.chmod('deploy_streamlit.sh', 0o755)
print("Streamlit app created successfully!")
    """
}

workflow STREAMLIT {
    take:
        zscore_file
        peptide_info
        pdb_dir
    
    main:
        GENERATE_JSMOL_HTML(zscore_file, peptide_info, pdb_dir)
        
        if (params.deploy_streamlit) {
            CREATE_STREAMLIT_APP(
                GENERATE_JSMOL_HTML.out.html_files,
                GENERATE_JSMOL_HTML.out.pdb_files,
                GENERATE_JSMOL_HTML.out.epitope_data
            )
        }
    
    emit:
        html_files = GENERATE_JSMOL_HTML.out.html_files
        pdb_files = GENERATE_JSMOL_HTML.out.pdb_files
        logs = GENERATE_JSMOL_HTML.out.logs
        streamlit_app = params.deploy_streamlit ? CREATE_STREAMLIT_APP.out.app : Channel.empty()
        deploy_script = params.deploy_streamlit ? CREATE_STREAMLIT_APP.out.deploy_script : Channel.empty()
}