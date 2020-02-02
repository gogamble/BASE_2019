#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
from functools import reduce
from collections import OrderedDict
import pandas as pd

## CHange qsim values to test
parser = argparse.ArgumentParser(description= "Change qsim")
parser.add_argument("-datafile","--datafile",dest="datafile",action="store",required=True,help="Provide path for bayesian input file")
parser.add_argument("-output","--output",dest="output",action="store",required=True,help="Provide path for bayesian input file with changed qsim")
args = parser.parse_args()

##Should add check if condiiton # variable doesn't match number specified in design file?

## Standardize Paths
args.output = os.path.abspath(args.output)
args.datafile = os.path.abspath(args.datafile)

# Read in design file and set to dataframe
df_data = pd.read_csv(args.datafile)

df_data.set_index('FEATURE_ID')

df_data['qsim_g1'].fillna(.49, inplace=True)
df_data['qsim_g2'].fillna(.49, inplace=True)
df_data['qsim_both'].fillna(.49, inplace=True)

df_data.loc[df_data.qsim_both == .49, 'c1_flag_analyze'] = 0
df_data.loc[df_data.qsim_both == .49, 'c2_flag_analyze'] = 0

bothlist = []
for header in df_data.columns[5:]:
    if 'both' not in header:
        bothlist.append(header)
print(bothlist)

both_df = df_data[bothlist]
zerolist = both_df.index[both_df.eq(0).all(1)].tolist()
print(zerolist)

for i in zerolist:
    data_df.loc[i, 'c1_flag_analyze'] = 0
    data_df.loc[i, 'c2_flag_analyze'] = 0

df_data.to_csv(args.output, index=False)
