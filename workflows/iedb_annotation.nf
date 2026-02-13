nextflow.enable.dsl = 2

params.iedb_database = "$baseDir/data/VIR3-example/epitope_full_v3.tsv"
params.epitope_length_min = 4    
params.epitope_length_max = 35
params.zscore_threshold = 3.5
params.max_workers = 4  
params.chunk_size = 100  

process EXTRACT_SAMPLE_NAMES {
    input:
    path(zscore_file)
    
    output:
    stdout emit: sample_names
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import gzip
    
    zscore_file = "${zscore_file}"
    
    if zscore_file.endswith('.gz'):
        df = pd.read_csv(gzip.open(zscore_file), index_col=0, nrows=1)
    else:
        df = pd.read_csv(zscore_file, index_col=0, nrows=1)
    
    for col in df.columns:
        print(col)
    """
}

process IEDB_ANNOTATION {

    publishDir "${params.results}/iedb_annotation", mode: 'copy'
    cpus params.max_workers
    memory '7 GB'
    
    input:
    tuple path(zscore_file), path(peptide_table), path(iedb_db), path(split_hits_file), val(sample_name)

    output:
    path "${sample_name}_annotated_peptides.csv", emit: annotated
    path "${sample_name}_annotation_summary.txt", emit: summary
    path "${sample_name}_novel_peptides.csv", emit: novel
    path "${sample_name}_significant_epitopes.csv", emit: significant

    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import numpy as np
    import time
    import re
    import gzip
    import gc
    
    start_time = time.time()
    sample_id = "${sample_name}"
        
    peptide_df = pd.read_csv('${peptide_table}')
    
    zscore_file = "${zscore_file}"
    
    if zscore_file.endswith('.gz'):
        zscore_df = pd.read_csv(gzip.open(zscore_file), index_col=0)
    else:
        zscore_df = pd.read_csv(zscore_file, index_col=0)
        
    if sample_id in zscore_df.columns:
        sample_zscore_df = zscore_df[[sample_id]].copy()
    else:
        sample_zscore_df = pd.DataFrame(index=zscore_df.index, columns=[sample_id])
        sample_zscore_df[sample_id] = 0
    
    split_hits = {}
    split_hits_file = "${split_hits_file}"
    try:
        if split_hits_file.endswith('.gz'):
            hits_df = pd.read_csv(gzip.open(split_hits_file))
        else:
            hits_df = pd.read_csv(split_hits_file)
        
        if 'id' in hits_df.columns:
            hit_columns = [col for col in hits_df.columns if col != 'id']
            if hit_columns:
                hit_col = hit_columns[0]  # Should be the sample column
                hits_df[hit_col] = hits_df[hit_col].fillna(0).astype(bool)
                split_hits = dict(zip(hits_df['id'], hits_df[hit_col]))
            else:
                print("No hit columns found in split hits file")
        else:
            print("No 'id' column found in split hits file")
    except Exception as e:
        split_hits = {}
    
    merged_df = peptide_df.merge(sample_zscore_df, left_index=True, right_index=True, how='inner')
    
    merged_df['zscore'] = merged_df[sample_id]
    merged_df['max_zscore'] = merged_df[sample_id]
    merged_df['mean_zscore'] = merged_df[sample_id]
    merged_df['significant_samples'] = (merged_df[sample_id] >= ${params.zscore_threshold}).astype(int)
    
    if split_hits:
        merged_df['has_virus_hit'] = merged_df.index.map(split_hits).fillna(False)
        virus_hit_count = merged_df['has_virus_hit'].sum()
    else:
        merged_df['has_virus_hit'] = False
    
    del zscore_df, sample_zscore_df, peptide_df
    gc.collect()
    
    chunk_structures = []
    try:
        chunk_reader = pd.read_csv('${iedb_db}', sep='\\t', chunksize=1000, low_memory=False)
        for i, chunk in enumerate(chunk_reader):
            if i >= 5:
                break
            chunk_structures.append({
                'chunk_num': i,
                'columns': list(chunk.columns),
                'num_cols': len(chunk.columns),
                'sample_data': chunk.iloc[0].tolist() if len(chunk) > 0 else []
            })
            if i == 0:
                print(f"      Columns: {list(chunk.columns)[:5]}...")
        
        del chunk_reader
    except Exception as e:
        chunk_structures = []
    
    best_col_idx = None
    if chunk_structures:
        epitope_column_candidates = {}
        
        for struct in chunk_structures:
            chunk_num = struct['chunk_num']
            columns = struct['columns']
            
            try:
                test_chunk = pd.read_csv('${iedb_db}', sep='\\t', 
                                       skiprows=chunk_num*1000, nrows=1000, 
                                       header=None, low_memory=False)
                
                for col_idx, col_name in enumerate(columns):
                    if col_idx < len(test_chunk.columns):
                        col_data = test_chunk.iloc[:, col_idx].dropna().astype(str)
                        
                        epitope_count = 0
                        for val in col_data.head(20):
                            val_clean = str(val).upper().strip()
                            if (${params.epitope_length_min} <= len(val_clean) <= ${params.epitope_length_max} and
                                len(re.findall(r'[ACDEFGHIKLMNPQRSTVWY]', val_clean)) >= len(val_clean) * 0.8 and
                                not any(pattern in val_clean.lower() for pattern in 
                                       ['http', 'linear peptide', 'object type', 'www', 'iedb'])):
                                epitope_count += 1
                        
                        if epitope_count >= 3:  
                            key = f"col_{col_idx}"
                            if key not in epitope_column_candidates:
                                epitope_column_candidates[key] = []
                            epitope_column_candidates[key].append({
                                'chunk': chunk_num,
                                'col_name': col_name,
                                'epitope_count': epitope_count
                            })
                            
            except Exception as e:
                continue
        
        best_score = 0
        for col_key, occurrences in epitope_column_candidates.items():
            score = sum(occ['epitope_count'] for occ in occurrences)
            if score > best_score:
                best_score = score
                best_col_idx = int(col_key.split('_')[1])
    
    epitope_set = set()
    total_processed = 0
    
    if best_col_idx is not None:
        print(f"Extracting epitopes from column position {best_col_idx}...")
        
        try:
            chunk_size = 5000            
            chunk_reader = pd.read_csv('${iedb_db}', sep='\\t', chunksize=chunk_size,
                                     header=None, low_memory=False)
            
            for chunk_num, chunk_df in enumerate(chunk_reader):
                total_processed += len(chunk_df)
                
                if best_col_idx >= len(chunk_df.columns):
                    continue
                
                values = chunk_df.iloc[:, best_col_idx].dropna().astype(str)
                
                for val in values:
                    val_clean = str(val).upper().strip()
                    
                    if (${params.epitope_length_min} <= len(val_clean) <= ${params.epitope_length_max} and
                        len(re.findall(r'[ACDEFGHIKLMNPQRSTVWY]', val_clean)) >= len(val_clean) * 0.8 and
                        not any(pattern in val_clean.lower() for pattern in 
                               ['http', 'linear peptide', 'object type', 'www', 'iedb'])):
                        epitope_set.add(val_clean)

        except Exception as e:
            print(f"ERROR reading IEDB database: {e}")
    
    epitope_list = list(epitope_set)
    
    epitope_to_organism = {}
    if best_col_idx is not None:
        try:
            chunk_reader = pd.read_csv('${iedb_db}', sep='\\t', chunksize=5000,
                                     header=None, low_memory=False)
            
            organism_col_idx = 13
            
            for chunk_num, chunk_df in enumerate(chunk_reader):
                if best_col_idx >= len(chunk_df.columns):
                    continue
                
                epitope_values = chunk_df.iloc[:, best_col_idx].dropna().astype(str)
                
                if organism_col_idx < len(chunk_df.columns):
                    organism_values = chunk_df.iloc[:, organism_col_idx].fillna('Unknown').astype(str)
                else:
                    organism_values = pd.Series(['Unknown'] * len(epitope_values))
                
                for epitope_val, org_val in zip(epitope_values, organism_values):
                    epitope_clean = str(epitope_val).upper().strip()
                    if (${params.epitope_length_min} <= len(epitope_clean) <= ${params.epitope_length_max} and
                        len(re.findall(r'[ACDEFGHIKLMNPQRSTVWY]', epitope_clean)) >= len(epitope_clean) * 0.8 and
                        not any(pattern in epitope_clean.lower() for pattern in 
                               ['http', 'linear peptide', 'object type', 'www', 'iedb'])):
                        if epitope_clean not in epitope_to_organism:
                            epitope_to_organism[epitope_clean] = str(org_val).strip()
                    
        except Exception as e:
            print(f"Warning: Could not extract organism info: {e}")
    
    def extract_sequence_function(row):
        if 'sequence' in row:
            return str(row['sequence']).upper().strip()
        elif 'peptide' in row:
            return str(row['peptide']).upper().strip()
        else:
            for col in row.index:
                if 'seq' in col.lower():
                    return str(row[col]).upper().strip()
            return ""

    def find_epitope_matches(peptide_seq, epitope_set, min_len=${params.epitope_length_min}, max_len=${params.epitope_length_max}):
        matches = []
        peptide_upper = peptide_seq.upper()
        
        for length in range(min_len, min(max_len + 1, len(peptide_seq) + 1)):
            for start in range(len(peptide_seq) - length + 1):
                substring = peptide_upper[start:start + length]
                if substring in epitope_set:
                    matches.append({
                        'epitope': substring,
                        'start': start,
                        'end': start + length
                    })
        
        return matches

    epitope_set_fast = set(epitope_list)
    all_annotations = []
    all_matches = []
    processed_count = 0

    merged_df['extracted_sequence'] = merged_df.apply(extract_sequence_function, axis=1)
    valid_peptides = merged_df[merged_df['extracted_sequence'].str.len() >= ${params.epitope_length_min}]

    for idx, row in valid_peptides.iterrows():
        peptide_seq = row['extracted_sequence']
        peptide_id = idx
        
        peptide_organism = row.get('Organism', 'Unknown')
        
        matches = find_epitope_matches(peptide_seq, epitope_set_fast)
        processed_count += 1
        
        annotation = {
            'sample_id': sample_id,
            'peptide_id': peptide_id,
            'peptide_sequence': peptide_seq,
            'organism': peptide_organism,
            'zscore': row.get('zscore', 0),
            'has_virus_hit': row.get('has_virus_hit', False),
            'num_epitope_matches': len(matches),
            'is_significant': row.get('zscore', 0) >= ${params.zscore_threshold},
            'is_novel': len(matches) == 0
        }
        all_annotations.append(annotation)
        
        for match in matches:
            is_significant = row.get('zscore', 0) >= ${params.zscore_threshold}
            epitope_seq = match['epitope']
            source_organism = epitope_to_organism.get(epitope_seq, 'Unknown')
            
            match_record = {
                'sample_id': sample_id,
                'peptide_id': peptide_id,
                'peptide_sequence': peptide_seq,
                'peptide_organism': peptide_organism,
                'peptide_zscore': row.get('zscore', 0),
                'peptide_has_virus_hit': row.get('has_virus_hit', False),
                'epitope_sequence': epitope_seq,
                'epitope_source_organism': source_organism,
                'start_position': match['start'],
                'end_position': match['end'],
                'is_significant': is_significant
            }
            all_matches.append(match_record)
    
    annotated_df = pd.DataFrame(all_annotations)
    matches_df = pd.DataFrame(all_matches)
    
    if not annotated_df.empty:
        novel_df = annotated_df[annotated_df['is_novel'] == True].copy()
        
        if not novel_df.empty:
            novel_df['potential_interest'] = (
                (novel_df['zscore'] >= ${params.zscore_threshold}) &
                (novel_df['has_virus_hit'] == True)
            )
        
        if not matches_df.empty:
            significant_df = matches_df[matches_df['is_significant'] == True].copy()
        else:
            significant_df = pd.DataFrame()
    else:
        novel_df = pd.DataFrame()
        significant_df = pd.DataFrame()
    
    annotated_df.to_csv(f'{sample_id}_annotated_peptides.csv', index=False)
    novel_df.to_csv(f'{sample_id}_novel_peptides.csv', index=False)
    significant_df.to_csv(f'{sample_id}_significant_epitopes.csv', index=False)
    
    total_time = time.time() - start_time
    matched_count = len(annotated_df[annotated_df['num_epitope_matches'] > 0]) if not annotated_df.empty else 0
    novel_count = len(novel_df)
    novel_significant = len(novel_df[novel_df['is_significant'] == True]) if not novel_df.empty else 0
    novel_with_virus_hit = len(novel_df[novel_df['has_virus_hit'] == True]) if not novel_df.empty else 0
    
    with open(f'{sample_id}_annotation_summary.txt', 'w') as f:
        f.write(f"IEDB Epitope Annotation Summary - Sample: {sample_id}\\n")
        f.write(f"=========================================================\\n\\n")
                
        f.write(f"Sample Results:\\n")
        f.write(f"  - Sample ID: {sample_id}\\n")
        f.write(f"  - Total peptides: {len(merged_df)}\\n")
        f.write(f"  - Peptides processed: {processed_count}\\n")
        f.write(f"  - Peptides with IEDB matches: {matched_count}\\n")
        f.write(f"  - Novel peptides: {novel_count}\\n")
        f.write(f"  - Novel significant peptides: {novel_significant}\\n")
        f.write(f"  - Novel peptides with virus hits: {novel_with_virus_hit}\\n")
        if processed_count > 0:
            f.write(f"  - IEDB match rate: {(matched_count/processed_count*100):.1f}%\\n")
            f.write(f"  - Novel rate: {(novel_count/processed_count*100):.1f}%\\n\\n")
        
        f.write(f"Significance Analysis:\\n")
        f.write(f"  - Significant matches: {len(significant_df)}\\n\\n")
        
        f.write(f"IEDB Database:\\n")
        f.write(f"  - Epitopes extracted: {len(epitope_list)}\\n")
        f.write(f"  - Database rows processed: {total_processed}\\n\\n")
        
        if not novel_df.empty and 'potential_interest' in novel_df.columns:
            high_interest = len(novel_df[novel_df['potential_interest'] == True])
            f.write(f"Novel Peptides of Interest:\\n")
            f.write(f"  - High interest (significant + virus hit): {high_interest}\\n")
    
    print(f"Novel peptides: {novel_count} (with virus hits: {novel_with_virus_hit})")
    """
}

workflow IEDB {
    take:
    zscore_files
    peptide_table

    main:
    iedb_db = file(params.iedb_database)
    
    consensus_zscore = zscore_files
        .flatten()
        .first()
        .map { file(it) }  
    
    EXTRACT_SAMPLE_NAMES(consensus_zscore)
    sample_names = EXTRACT_SAMPLE_NAMES.out.sample_names
        .splitText()
        .map { it.trim() }
        .filter { it != '' }
    
    split_hits_files = Channel.fromPath("${params.results}/virusscore/split_hits/hits/*.csv.gz")
        .toSortedList { a, b -> 
            def aNum = a.baseName.replaceAll(/\.csv/, '') as Integer
            def bNum = b.baseName.replaceAll(/\.csv/, '') as Integer
            return aNum <=> bNum 
        }
        .flatten()
    
    sample_input = sample_names
        .merge(split_hits_files)
        .combine(consensus_zscore)
        .map { sample_name, hits_file, zscore_file ->
            tuple(zscore_file, peptide_table, iedb_db, hits_file, sample_name)
        }
    
    IEDB_ANNOTATION(sample_input)

    emit:
    IEDB_ANNOTATION.out.annotated.collect()
}