#!/usr/bin/env python
import csv, sys, argparse, os, itertools, operator, collections
import pandas as pd
import numpy as np

## In alignment summary file, total_reads = sum(count_mapped_read_oposite_strand + count_unmapped_read + count_mapped_read + count_ambiguous_read)

parser = argparse.ArgumentParser(description='check total reads column in alignment file match the sum of the other columns')
parser.add_argument('-a1','--alnSum1',dest='alnSum1', action='store', required=True, help='The G1 alignment summary file containing all read types [Required]')
parser.add_argument('-a2','--alnSum2',dest='alnSum2', action='store', required=True, help='The G2 alignment summary file containing all read types [Required]')
parser.add_argument('-numread1','--numread1', dest='numread1', action='store', required=True, help='Input file with G1 pre-alignment read count data [Required]')
parser.add_argument('-numread2','--numread2', dest='numread2', action='store', required=True, help='Input file with G2 pre-alignment read count data [Required]')
parser.add_argument('-fq','--fq', dest='fq', action='store', required=True, help='FQ name [Required]')
parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file containing check info [Required]')
args = parser.parse_args()

## open G1 pre-alignment read count file
with open(args.numread1, 'r') as pre_aln_read1:
    read1=csv.reader(pre_aln_read1, delimiter=',')
    next(read1)
    for row in read1:
        start_reads1=int(row[1])

## open G2 pre-alignment read count file
with open(args.numread2, 'r') as pre_aln_read2:
    read2=csv.reader(pre_aln_read2, delimiter=',')
    next(read2)
    for row in read2:
        start_reads2=int(row[1])

## open G1 aln summary file and output file
with open(args.alnSum1, 'r') as sum_table1:
    sumT1=csv.reader(sum_table1, delimiter=',')
    next(sumT1)
    for row in sumT1:
        opp1=int(row[2])
        unmap1=int(row[3])
        mapread1=int(row[4])
        amb1=int(row[5])
        end_reads1 = opp1 + unmap1 + mapread1 + amb1

with open(args.alnSum2, 'r') as sum_table2:
    sumT2=csv.reader(sum_table2, delimiter=',')
    next(sumT2)
    for row in sumT2:
        opp2=int(row[2])
        unmap2=int(row[3])
        mapread2=int(row[4])
        amb2=int(row[5])
        end_reads2 = opp2 + unmap2 + mapread2 + amb2

if start_reads1 == end_reads1:
    flag_start_eq_end_G1 = 1
else:
    flag_start_eq_end_G1 = 0

if start_reads2 == end_reads2:
    flag_start_eq_end_G2 = 1
else:
    flag_start_eq_end_G2 = 0

with open(args.out, 'w') as outfile:
    spamwriter=csv.writer(outfile, delimiter=',')
    first_row = True
    
    if first_row:
        spamwriter.writerow(['fqName', 'start_reads1', 'start_reads2', 'end_reads1', 'end_reads2', 'flag_start_eq_end_G1', 'flag_start_eq_end_G2'])
        first_row=False

        row_items = [args.fq, start_reads1, start_reads2, end_reads1, end_reads2, flag_start_eq_end_G1, flag_start_eq_end_G2]
    spamwriter.writerow(row_items)




