#!/usr/bin/env python
from __future__ import division

import sys
import gzip
import difflib
import argparse
import pandas as pd

def is_novel_peptide(peptide, assigned_peptides, epitope_len):
    for assigned_peptide in assigned_peptides:        
        matcher = difflib.SequenceMatcher(None, peptide, assigned_peptide)
        match = matcher.find_longest_match(0, len(peptide), 0, len(assigned_peptide))
        if match.size >= epitope_len:
            return False
    return True

def calc_virus_scores(series, level, epitope_len):
    nhits_per_virus = series.groupby(level=level).sum()
    
    virus_scores = pd.Series(index=nhits_per_virus.index, name=series.name).fillna(0).astype(int)

    assigned_peptides = set()
    grouped = series[series].groupby(level=level)
    
    for virus in nhits_per_virus[nhits_per_virus > 0].sort_values(ascending=False).index:
        virus_hits = grouped.get_group(virus)

        score = 0
        peptides = virus_hits.index.get_level_values('peptide')
        for peptide in peptides:
            if is_novel_peptide(peptide, assigned_peptides, epitope_len):
                score += 1
                assigned_peptides.add(peptide)
        virus_scores[virus] = score
    return virus_scores

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate virus scores based on peptide hits.")
    parser.add_argument("hits", help="Path to the hits file (gzip compressed)")
    parser.add_argument("oligo_metadata", help="Path to the oligo metadata file (gzip compressed)")
    parser.add_argument("level", help="Level at which to group the data for scoring")
    parser.add_argument("epitope_len", type=int, help="Length of the epitope for novelty check")

    args = parser.parse_args()

    try:
        hits = pd.read_csv(gzip.open(args.hits), index_col=0).squeeze("columns")
        lib = pd.read_csv(args.oligo_metadata, low_memory=False)
    except Exception as e:
        sys.stderr.write(f"Error reading files: {e}\n")
        sys.exit(1)

    columns = ['id', 'Species', 'Organism', 'Entry', 'peptide']
    hits2 = pd.merge(lib[columns], hits.reset_index(), on='id').set_index(columns)

    virus_scores = calc_virus_scores(hits2[hits.name], args.level, args.epitope_len)
    virus_scores.to_csv(sys.stdout, header=True)