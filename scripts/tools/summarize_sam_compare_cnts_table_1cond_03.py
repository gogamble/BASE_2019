#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd
import csv 
import numpy as np
from functools import reduce

DEBUG = False

##  Prep sam compare data for Bayesian Machine
    ## imputs:
    ## design file containing sampleID,comparison_1 and comparison_2
    ## ase count tables from sam compare after matrix addition of tech reps and apn estimation
    ## user specified APN threshold for flagging feature as 'expressed' in reach rep 
    ## output:
    ## table containing 

    ## logic:
    ##  group df by sample, 
        ## only keep samples with more than 2 reps
    ##  summarize columns
    ##  filtering
            ## calculate:
        # total_reads_counted
        # both total
        # g1_total
        # g2 total
        # ase total
        ## for each rep, if APN > input then flag_APN = 1 (flag_APN = 0 if APN < input, flag_APN = -1 if APN < 0)
        ## if flag_APN = 1 for at least 1 of the reps then flag_analyze = 1  
    ##  merge reps together

def getOptions():
    parser = argparse.ArgumentParser(description='Return best row in blast scores file')
    parser.add_argument("-o", "--output", dest="output", action="store", required=True, help="Output directory for filtered ase counts")
    parser.add_argument("-d", "--design", dest="design", action='store', required=True, help="Design file")
    parser.add_argument("-c", "--compCol", dest="compCol", action='store', required=True, help="Name of the comparison column in your design file")
    parser.add_argument("-l", "--lineCol", dest="lineCol", action='store', required=True, help="Name of the column in your design file containing lines")
    parser.add_argument("-s", "--sam-compare-dir", dest="samCompDir", action='store', required=True, help="Path to directory containing summed ase count tables")
    parser.add_argument("-a", "--apn", dest="apn", action='store', required=True, type=int, help="APN (average per nucleotide) value for flagging a feature as found and analyzable")

    parser.add_argument("--debug", action='store_true', default=False, help="Print debugging output")


    args=parser.parse_args()
    return(args)

def main():
    """Main Function"""
    args = getOptions()
    global DEBUG
    if args.debug: DEBUG = True

    ## Read in design file as dataframe
    df_design = pd.read_csv(args.design, sep=',', header=0)
    # if DEBUG:

# set what column to use for comparison
    comparison_1 = args.compCol

# set what column to use for grouping lines
    line = args.lineCol

# groupby G2 and comparison    
    df_design['numReps'] = df_design.groupby([line, comparison_1]).Comparate1.transform('count')

#counter variables for each comparate 
    df_design['seqCnt'] = df_design.groupby(['G2', 'Comparate2']).cumcount()+1

# initialize dictionary to store count tables 
    df_dict={}
    df_dict2={}
    df_grouped=df_design.loc[df_design.numReps > 1]
#    print(df_grouped)
    F1IDlist=[]
    compcount=-1
    for index, sample in df_grouped.iterrows():
        #count_good variable stores rep number for each line-comparison
        compcount = compcount + 1
        g1 = sample['G1']
        g2 = sample['G2']
        count_good = sample['seqCnt']
        print(count_good)
        sample_id = sample['Comparate1']
        print('sample_id ' + sample_id)
        sample_id2 = sample_id.rsplit('_', 1)[0]
        print('sample_id2 ' + sample_id2)

        F1ID = g2 + '_' + g1 + '_sim'
        F1IDlist.append(F1ID)

        comp=sample_id2
        numReps = sample['numReps']
        repCnt = comp + "_num_reps"
        g1_total = comp + "_" + "g1_total" + "_" + "rep" + str(count_good)
        g2_total = comp + "_" + "g2_total" + "_" + "rep" + str(count_good)
        both_total = comp + "_" + "both_total" + "_" + "rep" + str(count_good)
        flag_apn = comp + "_" + "flag_apn" +"_" + "rep" + str(count_good)

# create key to store count table in dictionary
        countKey  = sample_id 

        print('countKey ' + countKey)

# read in count table     
        inFile = 'ase_sum_counts_' + sample_id + '.csv'
        print(inFile)
        samC = os.path.join(args.samCompDir, inFile)
        count_missing_sample_file = 0
        try:
            samFile = pd.read_csv(samC, sep=',', header=0) 
        except:
            print(f"Missing:\n {samC}")
            count_missing_sample_file += 1
            continue

        ## summarize
        print('samC' + samC)
        samFile[both_total] = samFile['BOTH_EXACT'] + samFile['BOTH_INEXACT_EQUAL']
        samFile[g2_total] = samFile['SAM_A_ONLY_EXACT'] + samFile['SAM_A_ONLY_SINGLE_INEXACT'] + samFile['SAM_A_EXACT_SAM_B_INEXACT'] + samFile['SAM_A_INEXACT_BETTER']
 
        samFile[g1_total] = samFile['SAM_B_ONLY_EXACT'] + samFile['SAM_B_ONLY_SINGLE_INEXACT'] + samFile['SAM_B_EXACT_SAM_A_INEXACT'] + samFile['SAM_B_INEXACT_BETTER']

        ## if APN > user-specified value then flag_APN = 1
        samFile[flag_apn] = 0
        samFile.loc[samFile.APN_total_reads > args.apn, flag_apn] = 1
        samFile.loc[samFile.APN_total_reads  == 0 , flag_apn] = 0
        samFile.loc[samFile.APN_total_reads < 0 , flag_apn] = -1
        samFile[flag_apn] = pd.to_numeric(samFile[flag_apn], errors='ignore')                

        sam_subset = samFile[['FEATURE_ID', g1_total, g2_total, both_total, flag_apn]]
        sam_subset['g1'] = g1
        sam_subset['g2'] = g2
        sam_subset[repCnt] = numReps

        df_dict[countKey] = sam_subset

        comp_list = [value for key,value in df_dict.items() if key.startswith(comp)]

# need these column headers to merge on
        comp_reps = comp + "_num_reps"

# merge c1 dataframes and c2 dataframes separately
        comp_merged = reduce(lambda x,y: pd.merge(x,y, on= ['FEATURE_ID','g1', 'g2', comp_reps], how='outer'), comp_list)

# create flag_analyze for each comparate, pull out all flag_apn using pattern match
#    sum flag_apn for each comparate, if flag_sum  > 0, set new column comparison_flag_analyze  = 1

        comp_flag_analyze = comp + "_flag_analyze"

        flag_comp_cols = [col for col in comp_merged.columns if 'flag_apn' in col]
        comp_merged['flag_sum'] = comp_merged[flag_comp_cols].sum(axis=1)
        comp_merged[comp_flag_analyze] = np.where(comp_merged['flag_sum'] > 0, 1, 0)
        del comp_merged['flag_sum']

# order column headers for output
        headers = list(comp_merged.columns.values)
        comp_reps_headers = []

        for each in headers:
            if comp in each and 'rep' in each and 'num' not in each:
                comp_reps_headers.append(each)
        ordered_headers = ['FEATURE_ID', 'g1', 'g2', comp+'_flag_analyze', comp+'_num_reps'] + comp_reps_headers

        outfilename = "ase_counts_filtered_" + comp + ".csv"
        outfile = os.path.join(args.output, outfilename)

        comp_merged = comp_merged[ordered_headers]
        comp_merged.to_csv(outfile, index=False)
        print(outfile)
    print(F1IDlist)
    df_drop = df_design.drop(['numReps','seqCnt'], axis=1)
    df_drop['F1ID']=F1IDlist
    df_drop.to_csv(args.design, index=False)

if __name__=='__main__':
    main()
