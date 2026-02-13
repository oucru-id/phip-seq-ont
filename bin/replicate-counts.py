#!/usr/bin/env python

import pandas as pd
import numpy as np
import phippery
from phippery.utils import load, dump, get_annotation_table
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-ds", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()


def replicate_oligo_counts(ds, peptide_oligo_feature="Oligo"):

    pep_anno_table = get_annotation_table(ds, "peptide")

    for oligo_seq, pep_anno_table_oligo in pep_anno_table.groupby(peptide_oligo_feature):

        if pep_anno_table_oligo.shape[0] == 1:
            continue
        idxs = pep_anno_table_oligo.index.values
        ds.counts.loc[idxs, :] = np.tile(
                ds.counts.loc[idxs, :].sum(axis=0).values,
                (pep_anno_table_oligo.shape[0], 1)
        )

ds = phippery.load(args.ds)
replicate_oligo_counts(ds, "oligo")
phippery.dump(ds, args.o)