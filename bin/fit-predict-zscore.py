#!/usr/bin/env python

from phippery.utils import * 
from phippery.modeling import zscore

import argparse
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("-ds", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()

ds = load(args.ds)
beads_ds = ds_query(ds, "control_status == 'beads_only'")

zscore_ds = zscore(
    ds,
    beads_ds,
    data_table='cpm',
    min_Npeptides_per_bin=300,
    lower_quantile_limit=0.05,
    upper_quantile_limit=0.95,
    inplace=False,
    new_table_name='zscore'
)

dump(zscore_ds, args.o)