#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
from functools import reduce
from collections import OrderedDict
import pandas as pd

### Comparate names are kept in individual comparate file headers
### Make more general (i.e. c1, c2, c3...) so its easier

parser = argparse.ArgumentParser(description= "Change comparate header names to c1, c2,...")
parser.add_argument("-design","--design",dest="design",action="store",required=True,help="Design file containing G1, G2, comparate names, and compIDs for files to analyze[CSV]")
parser.add_argument("-datafile","--datafile",dest="datafile",action="store",required=True,help="Provide path for input datafile [CSV]")
parser.add_argument("-array","--array",dest="array",action="store", required=True,help="Numerical value representing row of design file to be analyzed")
parser.add_argument("-cond","--cond",dest="cond",action="store", required=True,help="Numerical value representing number of comparates to be analyzed")
parser.add_argument("-o","--output",dest="output",action="store", required=True,help="Output for final bayesian input file[CSV]")
args = parser.parse_args()

##Should add check if condiiton # variable doesn't match number specified in design file?

## Standardize Paths 
args.output = os.path.abspath(args.output)
args.design = os.path.abspath(args.design)

# Read in design file and set to dataframe
df_design = pd.read_csv(args.design)

## index compates out of a single row of design file to analyze
Row=int(args.array) - 2
print(Row)

print(df_design)

df_design.set_index('G1')
data_row=df_design.iloc[Row]

print('data_row1')
print(data_row)

del df_design['G1']
del df_design['G2']

comparison=data_row['compID']

del data_row['compID']
del data_row['G1']
del data_row['G2']

data_row = data_row.to_frame()

if 'F1ID' in data_row.columns:
    del data_row['F1ID']

print('data_row2')
print(data_row)

row_list = data_row.values.tolist()

infileName = args.datafile + "bayesian_input_comp_" + comparison + ".csv"
print(infileName)

if os.path.isfile(infileName):

    data_df = pd.read_csv(infileName)
    data_df.set_index('FEATURE_ID')

    pre_headers=list(data_df.columns.get_values())
    pre_headers_split=pre_headers[:3]
    pre_headers_split2=pre_headers[3:]

    for index,val in enumerate(row_list):
        c = 'c' + str(index + 1)
        val = ''.join(val)

        for i in range(len(pre_headers_split2)):
            if val in pre_headers_split2[i]:
                pre_headers_split2[i] = pre_headers_split2[i].replace(val, c)
    pre_headers_cat = pre_headers_split + pre_headers_split2
    print(pre_headers_cat)

    for i in range(len(pre_headers_cat)):
        if 'comparison' in pre_headers_cat[i]:
            pre_headers_cat[i] = pre_headers_cat[i].replace('comparison', 'line')

    data_df.columns=pre_headers_cat
    r_input = args.output + "/bayesian_input_comp_" + comparison + '_temp.csv'

    print(r_input)

    data_df.to_csv(r_input, na_rep = 'NA', index=False)


