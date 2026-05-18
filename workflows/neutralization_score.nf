nextflow.enable.dsl = 2

process NEUTRALIZATION_SCORE {
    publishDir "${params.results}/neutralization_prediction", mode: 'copy'
    cpus params.max_workers

    input:
    path iedb_significant_epitopes
    path iedb_annotated_peptides
    path peptide_table
    path zscore_files
    path pdb_dir
    path neutralization_db

    output:
    path "neutralization_scores_per_sample.csv",   emit: scores
    path "neutralization_summary.txt",             emit: summary
    path "high_confidence_candidates.csv",         emit: candidates
    path "detailed_analysis.json",                 emit: detailed
    path "conformational_epitope_clusters.csv",    emit: clusters
    
    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import numpy as np
    import json
    import os
    import math
    import warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)

    HAS_BIOPYTHON = False
    try:
        from Bio.PDB import PDBParser, MMCIFParser
        from Bio.PDB.SASA import ShrakeRupley
        try:
            from Bio.Data.IUPACData import protein_letters_3to1
        except ImportError:
            from Bio.PDB.Polypeptide import protein_letters_3to1
        HAS_BIOPYTHON = True
    except ImportError:
        print("Warning: Biopython not found. Structural features will be skipped.")

    HAS_SCIPY = False
    try:
        from scipy.cluster.hierarchy import fclusterdata
        HAS_SCIPY = True
    except ImportError:
        print("Warning: SciPy not found. Spatial clustering will be skipped.")

    def read_csv_list(path_str):
        paths = path_str.strip().split()
        frames = []
        for p in paths:
            p = p.strip()
            if p and os.path.isfile(p):
                try:
                    frames.append(pd.read_csv(p))
                except Exception as e:
                    print(f"  Warning: could not read {p}: {e}")
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    print("Loading significant epitopes...")
    sig_epitopes = read_csv_list("${iedb_significant_epitopes}")

    print("Loading annotated peptides...")
    annot_peptides = read_csv_list("${iedb_annotated_peptides}")

    print("Loading peptide table...")
    peptide_df = pd.read_csv("${peptide_table}")

    print("Loading z-score matrix...")
    zscore_paths = "${zscore_files}".strip().split()
    zscore_file = next((f for f in zscore_paths if 'zscore' in f), None)
    if not zscore_file and zscore_paths:
        zscore_file = zscore_paths[0]
    if zscore_file and os.path.isfile(zscore_file):
        if zscore_file.endswith('.gz'):
            import gzip
            zscore_df = pd.read_csv(gzip.open(zscore_file), index_col=0)
        else:
            zscore_df = pd.read_csv(zscore_file, index_col=0)
        if zscore_df.empty:
            raise ValueError(f"ERROR: Z-score file loaded but is empty: {zscore_file}")
    else:
        raise FileNotFoundError(f"ERROR: No z-score file found. Searched paths: {zscore_paths}")

    VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')
    INVALID_PATTERNS = ['http', 'linear peptide', 'object type', 'www', 'iedb',
                        'discontinuous']

    def load_neutralization_db(filepath):
        seqs = set()
        if not os.path.isfile(filepath):
            print(f"  Warning: neutralization DB not found: {filepath}")
            return seqs
        try:
            df = pd.read_csv(filepath, sep='\\t', skiprows=2, header=None,
                             low_memory=False, on_bad_lines='skip')
            if df.shape[1] < 3:
                print(f"  Warning: unexpected format in {filepath}")
                return seqs
            col = df.iloc[:, 2].dropna().astype(str)
            for val in col:
                val = val.strip().upper()
                if not (4 <= len(val) <= 50):
                    continue
                aa_frac = sum(1 for c in val if c in VALID_AA) / max(len(val), 1)
                if aa_frac < 0.8:
                    continue
                vl = val.lower()
                if any(p in vl for p in INVALID_PATTERNS):
                    continue
                seqs.add(val)
        except Exception as e:
            print(f"  Warning: error reading {filepath}: {e}")
        return seqs

    print("Loading neutralization database...")
    neut_db = load_neutralization_db("${neutralization_db}")
    print(f"  Neutralization DB: {len(neut_db)} entries")

    MAX_SASA = {
        'A': 129, 'R': 274, 'N': 195, 'D': 193, 'C': 167,
        'E': 223, 'Q': 225, 'G': 104, 'H': 224, 'I': 197,
        'L': 201, 'K': 236, 'M': 224, 'F': 240, 'P': 159,
        'S': 155, 'T': 172, 'W': 285, 'Y': 263, 'V': 174,
    }

    def parse_structure_file(fpath, pdb_key):
        result = {}
        try:
            if fpath.endswith('.cif'):
                parser = MMCIFParser(QUIET=True)
            else:
                parser = PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_key, fpath)
            sr = ShrakeRupley()
            sr.compute(structure, level="R")
            model = next(iter(structure))
            for chain in model:
                residues = []
                raw_bfac = []
                for res in chain:
                    if res.get_id()[0] != ' ':   # skip HETATM
                        continue
                    try:
                        resname = res.get_resname().strip()
                        aa = protein_letters_3to1.get(resname.capitalize(),
                             protein_letters_3to1.get(resname, None))
                        if aa is None or aa not in MAX_SASA:
                            continue
                        ca = res['CA'] if 'CA' in res else None
                        if ca is None:
                            continue
                        sasa_norm = min(1.0, res.sasa / MAX_SASA[aa])
                        bfac_raw  = ca.get_bfactor()
                        ca_coord  = list(ca.get_vector())
                        residues.append({'aa': aa, 'sasa': sasa_norm,
                                         'bfac_raw': bfac_raw, 'ca': ca_coord})
                        raw_bfac.append(bfac_raw)
                    except (KeyError, Exception):
                        continue
                if not residues:
                    continue
                b_min, b_max = min(raw_bfac), max(raw_bfac)
                b_range = b_max - b_min if b_max > b_min else 1.0
                for r in residues:
                    r['bfac'] = (r['bfac_raw'] - b_min) / b_range
                result[chain.get_id()] = residues
        except Exception as e:
            print(f"  Warning: could not parse {fpath}: {e}")
        return result

    pdb_cache = {}   # {normalized_key: {chain_id: [res_dicts]}}
    pdb_dir_path = "${pdb_dir}"
    if HAS_BIOPYTHON and os.path.isdir(pdb_dir_path):
        print(f"Parsing PDB/CIF structures in {pdb_dir_path}...")
        for fname in os.listdir(pdb_dir_path):
            if fname.lower().endswith(('.pdb', '.cif')):
                fpath = os.path.join(pdb_dir_path, fname)
                key   = fname.lower().rsplit('.', 1)[0]
                chains = parse_structure_file(fpath, key)
                if chains:
                    chain_seqs = {cid: ''.join(r['aa'] for r in rlist)
                                  for cid, rlist in chains.items()}
                    pdb_cache[key] = {'chains': chains, 'chain_seqs': chain_seqs}
                    print(f"  Loaded {fname}: {len(chains)} chain(s)")
    else:
        if not HAS_BIOPYTHON:
            print("Skipping PDB parsing.")
        elif not os.path.isdir(pdb_dir_path):
            print(f"Skipping PDB parsing.")

    def get_structural_features(pep_seq, pdb_id):
        if not pdb_id or (isinstance(pdb_id, float) and math.isnan(pdb_id)):
            return None, None, None
        key = str(pdb_id).strip().lower().rsplit('.', 1)[0]
        if key not in pdb_cache:
            return None, None, None
        entry      = pdb_cache[key]
        chain_seqs = entry['chain_seqs']
        chains     = entry['chains']
        pep_up     = pep_seq.upper()
        anchor     = pep_up[:15] if len(pep_up) >= 15 else pep_up
        for cid, cseq in chain_seqs.items():
            idx = cseq.find(anchor)
            if idx == -1:
                continue
            rlist = chains[cid]
            end_idx = min(idx + len(pep_up), len(rlist))
            matched = rlist[idx:end_idx]
            if not matched:
                continue
            sasa_vals = [r['sasa'] for r in matched]
            bfac_vals = [r['bfac'] for r in matched]
            coords    = [r['ca']   for r in matched]
            centroid  = [sum(c[i] for c in coords) / len(coords) for i in range(3)]
            return (float(np.mean(sasa_vals)),
                    float(np.mean(bfac_vals)),
                    centroid)
        return None, None, None

    def build_sig_index(df):
        idx = {}
        if df.empty:
            return idx
        for _, row in df.iterrows():
            pid = row.get('peptide_id', None)
            if pid is None:
                continue
            if pid not in idx:
                idx[pid] = {'orgs': {}, 'positions': set()}
            org = str(row.get('epitope_source_organism', '')).strip()
            if org and org.lower() not in ('', 'nan', 'unknown'):
                idx[pid]['orgs'][org] = idx[pid]['orgs'].get(org, 0) + 1
            start = row.get('start_position', None)
            end   = row.get('end_position', None)
            if pd.notna(start) and pd.notna(end):
                try:
                    idx[pid]['positions'].update(range(int(start), int(end)))
                except (TypeError, ValueError):
                    pass
        return idx

    def shannon_entropy(org_counts):
        total = sum(org_counts.values())
        if total == 0:
            return 0.0
        entropy = 0.0
        for count in org_counts.values():
            p = count / total
            if p > 0:
                entropy -= p * math.log2(p)
        return entropy

    # Kyte-Doolittle hydropathy scale
    KD_SCALE = {
        'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
        'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
    }
    LOOP_RESIDUES = set('GPST')
    import re as _re
    _NGLYCOSYLATION_RE = _re.compile(r'N[^P][ST]')

    def build_annot_index(df):
        if df.empty or 'peptide_id' not in df.columns:
            return {}
        col = 'num_epitope_matches' if 'num_epitope_matches' in df.columns else None
        if col is None:
            return {}
        return df.groupby('peptide_id')[col].max().to_dict()

    def neutralization_db_score(pep_up):
        if pep_up in neut_db:
            return 3.0
        best = 0.0
        for ep in neut_db:
            if ep in pep_up:
                best = max(best, 1.5)
                continue
            if len(pep_up) >= 5 and len(ep) >= 5:
                pep_kmers = {pep_up[i:i+5] for i in range(len(pep_up) - 4)}
                ep_kmers  = {ep[i:i+5]     for i in range(len(ep) - 4)}
                if pep_kmers & ep_kmers:
                    best = max(best, 0.5)
        return best

    print("Building IEDB indices...")
    sig_idx   = build_sig_index(sig_epitopes)
    annot_idx = build_annot_index(annot_peptides)

    print("Pre-calculating static scores for all peptides...")
    static_scores = {}

    for _, row in peptide_df.iterrows():
        pid     = row.get('id', row.name)
        pep_seq = str(row.get('peptide', '')).upper().strip()
        org     = str(row.get('Organism', row.get('organism', ''))).strip()
        prot    = str(row.get('Protein',  row.get('protein',  ''))).lower().strip()
        pdb_id  = row.get('pdb_id', None)
        if pd.isna(pdb_id):
            pdb_id = None

        n_matches    = annot_idx.get(pid, 0)
        conservation = min(1.0 + math.log1p(n_matches) * 0.4, 3.0)

        org_counts   = sig_idx.get(pid, {}).get('orgs', {})
        n_orgs       = len(org_counts)
        entropy      = shannon_entropy(org_counts) if org_counts else 0.0
        iedb_evid    = min(1.0 + entropy * 0.6, 3.0)

        positions    = sig_idx.get(pid, {}).get('positions', set())
        pep_len      = max(len(pep_seq), 1)
        coverage     = min(len(positions) / pep_len, 1.0) * 3.0

        neut_score   = neutralization_db_score(pep_seq) if pep_seq else 0.0

        context_factor = 1.0
        if 'spike' in prot:
            context_factor = 1.12  
        elif 'envelope' in prot or 'membrane' in prot:
            context_factor = 1.06  

        gravy_vals = [KD_SCALE[aa] for aa in pep_seq if aa in KD_SCALE]
        gravy      = sum(gravy_vals) / max(len(gravy_vals), 1)
        gravy_mult = max(0.6, 1.0 - max(0.0, gravy) * 0.08)

        loop_frac   = sum(1 for aa in pep_seq if aa in LOOP_RESIDUES) / max(pep_len, 1)
        flex_bonus  = loop_frac * 0.3

        n_glyc_sites = len(_NGLYCOSYLATION_RE.findall(pep_seq))
        glycan_mult  = max(0.5, 1.0 - n_glyc_sites * 0.25)

        sasa_mean, bfac_mean, ca_centroid = get_structural_features(pep_seq, pdb_id)

        static_scores[pid] = {
            'seq':            pep_seq,
            'org':            org,
            'prot':           prot,
            'pdb_id':         pdb_id,
            'conservation':   conservation,
            'iedb_evid':      iedb_evid,
            'coverage':       coverage,
            'neut_score':     neut_score,
            'context_factor': context_factor,
            'gravy':          gravy,
            'gravy_mult':     gravy_mult,
            'flex_bonus':     flex_bonus,
            'glycan_mult':    glycan_mult,
            'n_glyc_sites':   n_glyc_sites,
            'sasa_mean':      sasa_mean,
            'bfac_mean':      bfac_mean,
            'ca_centroid':    ca_centroid,
        }

    print(f"Scoring across {len(zscore_df.columns)} samples...")

    W_CONS = ${params.conservation_weight}
    W_PHIP = ${params.phip_signal_weight}
    W_EVID = ${params.iedb_evidence_weight}
    W_NEUT = ${params.neutralization_db_weight}
    W_COV  = ${params.epitope_coverage_weight}
    W_BFAC = ${params.bfactor_weight}
    Z_THRESH = ${params.zscore_threshold}
    N_THRESH  = ${params.neutralization_threshold}

    results = []

    for sample_name in zscore_df.columns:
        sample_z = zscore_df[sample_name]
        hits = sample_z[sample_z >= Z_THRESH]

        for pid, zscore in hits.items():
            if pid not in static_scores:
                continue
            info = static_scores[pid]

            phip_q = 7.0 * (1.0 - math.exp(-zscore / 10.0))

            has_struct = (info['sasa_mean'] is not None and
                          info['bfac_mean'] is not None)

            if has_struct:
                bfac_comp = info['bfac_mean']
                composite = (
                    info['conservation'] * W_CONS +
                    phip_q               * W_PHIP +
                    info['iedb_evid']    * W_EVID +
                    info['neut_score']   * W_NEUT +
                    info['coverage']     * W_COV  +
                    bfac_comp            * W_BFAC
                )

                if info['sasa_mean'] < 0.05:
                    sasa_mult = 0.0
                else:
                    sasa_mult = max(0.1, min(1.0, 1.5 * info['sasa_mean']))
                composite *= sasa_mult
            else:
                denom = 1.0 - W_BFAC
                composite = (
                    info['conservation'] * (W_CONS / denom) +
                    phip_q               * (W_PHIP / denom) +
                    info['iedb_evid']    * (W_EVID / denom) +
                    info['neut_score']   * (W_NEUT / denom) +
                    info['coverage']     * (W_COV  / denom)
                )
                bfac_comp = None

            composite += info['flex_bonus']       
            composite *= info['gravy_mult']        
            composite *= info['glycan_mult']         
            composite *= info['context_factor']    

            if composite >= N_THRESH:
                prediction = 'High'
            elif composite >= N_THRESH * 0.65:
                prediction = 'Moderate'
            else:
                prediction = 'Low'

            results.append({
                'sample':                  sample_name,
                'peptide_id':              pid,
                'organism':                info['org'],
                'protein':                 info['prot'],
                'pdb_id':                  info['pdb_id'],
                'zscore':                  round(zscore, 3),
                'phip_signal':             round(phip_q, 3),
                'conservation_score':      round(info['conservation'], 3),
                'iedb_evidence_score':     round(info['iedb_evid'], 3),
                'neutralization_db_score': round(info['neut_score'], 3),
                'coverage_score':          round(info['coverage'], 3),
                'bfactor_score':           round(bfac_comp, 3) if bfac_comp is not None else None,
                'sasa_mean':               round(info['sasa_mean'], 3) if info['sasa_mean'] is not None else None,
                'gravy_score':             round(info['gravy'], 3),
                'flex_bonus':              round(info['flex_bonus'], 3),
                'n_glycosylation_sites':   info['n_glyc_sites'],
                'glycan_multiplier':       round(info['glycan_mult'], 3),
                'context_factor':          round(info['context_factor'], 3),
                'composite_score':         round(composite, 3),
                'prediction':              prediction,
                'has_structural_data':     has_struct,
                'cluster_id':              None,
            })

    cluster_rows = []

    if HAS_SCIPY and results:
        from collections import defaultdict
        groups = defaultdict(list)
        for r in results:
            if r['has_structural_data'] and r['pdb_id']:
                key = (r['sample'], str(r['pdb_id']))
                pid = r['peptide_id']
                coord = static_scores[pid]['ca_centroid']
                if coord:
                    groups[key].append((pid, coord))

        for (sample, pdb_key), members in groups.items():
            if len(members) < 2:
                continue
            pids   = [m[0] for m in members]
            coords = np.array([m[1] for m in members])
            try:
                labels = fclusterdata(coords,
                                      t=${params.spatial_cluster_angstrom},
                                      criterion='distance',
                                      metric='euclidean',
                                      method='single')
            except Exception as e:
                print(f"  Warning: clustering failed for ({sample},{pdb_key}): {e}")
                continue

            pid_to_label = dict(zip(pids, labels))
            for r in results:
                if r['sample'] == sample and str(r.get('pdb_id', '')) == pdb_key:
                    if r['peptide_id'] in pid_to_label:
                        r['cluster_id'] = f"{pdb_key}_c{pid_to_label[r['peptide_id']]}"

            from collections import Counter
            label_counts = Counter(labels)
            for label, count in label_counts.items():
                if count < 2:
                    continue
                cluster_pids    = [pids[i] for i, l in enumerate(labels) if l == label]
                cluster_coords  = coords[[i for i, l in enumerate(labels) if l == label]]
                cluster_center  = cluster_coords.mean(axis=0).tolist()
                cluster_rows.append({
                    'sample':       sample,
                    'pdb_id':       pdb_key,
                    'cluster_id':   f"{pdb_key}_c{label}",
                    'num_peptides': count,
                    'peptide_ids':  '|'.join(str(p) for p in cluster_pids),
                    'center_x':     round(cluster_center[0], 2),
                    'center_y':     round(cluster_center[1], 2),
                    'center_z':     round(cluster_center[2], 2),
                })

    results_df = pd.DataFrame(results)

    EMPTY_SCORE_COLS = ['sample', 'peptide_id', 'organism', 'protein', 'pdb_id',
                        'zscore', 'phip_signal', 'conservation_score',
                        'iedb_evidence_score', 'neutralization_db_score',
                        'coverage_score', 'bfactor_score', 'sasa_mean',
                        'gravy_score', 'flex_bonus', 'n_glycosylation_sites',
                        'glycan_multiplier', 'context_factor',
                        'composite_score', 'prediction', 'has_structural_data',
                        'cluster_id']

    if not results_df.empty:
        results_df = results_df.sort_values(
            ['sample', 'composite_score'], ascending=[True, False])
        results_df.to_csv('neutralization_scores_per_sample.csv', index=False)

        high_df = results_df[results_df['prediction'] == 'High'].copy()
        high_df.to_csv('high_confidence_candidates.csv', index=False)

        summary_counts = results_df['prediction'].value_counts().to_dict()
        struct_count   = int(results_df['has_structural_data'].sum())
        cluster_count  = int(results_df['cluster_id'].notna().sum())

        score_desc = results_df['composite_score'].describe()
        score_stats = {k: round(float(v), 3) for k, v in score_desc.items()}
        score_pcts  = {f'p{p}': round(float(results_df['composite_score'].quantile(p/100)), 3)
                       for p in [5, 10, 25, 50, 75, 90, 95]}

        with open('detailed_analysis.json', 'w') as f:
            json.dump({
                'total_hits_analyzed':     len(results_df),
                'samples_analyzed':        len(zscore_df.columns),
                'prediction_counts':       summary_counts,
                'peptides_with_structure': struct_count,
                'peptides_in_clusters':    cluster_count,
                'neutralization_db_size':  len(neut_db),
                'score_distribution':      score_stats,
                'score_percentiles':       score_pcts,
                'weights': {
                    'conservation':      ${params.conservation_weight},
                    'phip_signal':       ${params.phip_signal_weight},
                    'iedb_evidence':     ${params.iedb_evidence_weight},
                    'neutralization_db': ${params.neutralization_db_weight},
                    'epitope_coverage':  ${params.epitope_coverage_weight},
                    'bfactor':           ${params.bfactor_weight},
                },
                'algorithm_version': '2.0',
                'improvements': [
                    'bfactor_scaling_fixed',
                    'neut_db_score_logic_corrected',
                    'short_peptide_db_lookup_fixed',
                    'context_as_multiplicative_factor',
                    'zscore_empty_file_error_raised',
                    'sasa_hard_filter_buried_lt_0.05',
                    'n_glycosylation_penalty',
                    'gravy_hydropathy_multiplier',
                    'flexibility_loop_bonus',
                    'shannon_entropy_iedb_diversity',
                    'score_distribution_in_output',
                ],
            }, f, indent=2)

        with open('neutralization_summary.txt', 'w') as f:
            f.write("Neutralization Prediction\\n")
            f.write("===================================================\\n")
            f.write(f"Samples analyzed:            {len(zscore_df.columns)}\\n")
            f.write(f"Reactive peptides scored:    {len(results_df)}\\n")
            f.write(f"High-confidence candidates:  {len(high_df)}\\n")
            f.write(f"Moderate candidates:         {summary_counts.get('Moderate', 0)}\\n")
            f.write(f"Peptides with structure:     {struct_count}\\n")
            f.write(f"Neutralization DB entries:   {len(neut_db)}\\n")
            f.write(f"Peptides in spatial clusters:{cluster_count}\\n")
            f.write(f"Threshold used:              {N_THRESH}\\n")
            f.write(f"Score mean:                  {score_stats.get('mean', 'N/A')}\\n")
            f.write(f"Score std:                   {score_stats.get('std', 'N/A')}\\n")
            f.write(f"Score p50 (median):          {score_pcts.get('p50', 'N/A')}\\n")
            f.write(f"Score p95:                   {score_pcts.get('p95', 'N/A')}\\n")
    else:
        print("No significant hits found in any sample.")
        pd.DataFrame(columns=EMPTY_SCORE_COLS).to_csv(
            'neutralization_scores_per_sample.csv', index=False)
        pd.DataFrame(columns=EMPTY_SCORE_COLS).to_csv(
            'high_confidence_candidates.csv', index=False)
        with open('neutralization_summary.txt', 'w') as f:
            f.write("No hits found above z-score threshold.\\n")
        with open('detailed_analysis.json', 'w') as f:
            json.dump({'total_hits_analyzed': 0}, f)

    cluster_df = pd.DataFrame(cluster_rows) if cluster_rows else pd.DataFrame(
        columns=['sample', 'pdb_id', 'cluster_id', 'num_peptides',
                 'peptide_ids', 'center_x', 'center_y', 'center_z'])
    cluster_df.to_csv('conformational_epitope_clusters.csv', index=False)
    print(f"Conformational clusters found: {len(cluster_df)}")
    print("Done.")
    """
}

workflow NEUTRALIZATION_PREDICTION {
    take:
    iedb_significant_epitopes
    iedb_annotated_peptides
    peptide_table
    zscore_files
    pdb_dir
    neutralization_db

    main:
    NEUTRALIZATION_SCORE(
        iedb_significant_epitopes,
        iedb_annotated_peptides,
        peptide_table,
        zscore_files,
        pdb_dir,
        neutralization_db
    )

    emit:
    scores     = NEUTRALIZATION_SCORE.out.scores
    summary    = NEUTRALIZATION_SCORE.out.summary
    candidates = NEUTRALIZATION_SCORE.out.candidates
    detailed   = NEUTRALIZATION_SCORE.out.detailed
    clusters   = NEUTRALIZATION_SCORE.out.clusters
}
