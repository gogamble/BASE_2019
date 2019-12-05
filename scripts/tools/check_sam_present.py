#!/usr/bin/env python
import csv, sys, argparse, os, itertools, operator, collections
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='check read numbers into and out of sam compare')
parser.add_argument('-fq','--fq',dest='fq', action='store', required=True, help='Name of the fq file]')
parser.add_argument('-alnType','--alnType',dest='alnType', action='store', required=True, help='Input SE for single end or PE for paired end alignments [Required]')
parser.add_argument('-samPath','--samPath',dest='samPath', action='store', required=True, help='Path to sam files [Required]')
parser.add_argument('-G1','--G1',dest='G1', action='store', required=True, help='Tester [Required]')
parser.add_argument('-G2','--G2',dest='G2', action='store', required=True, help='Line [Required]')
parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file containing check info [Required]')
args = parser.parse_args()


### For every FQ file run, should have 2 sam files
### Append sam files to an array

fq_file = os.path.split(args.fq)[1]

samPath= args.samPath
sam1=args.samPath + '/' + args.G2 + '_' + args.fq + '_upd_feature_uniq.sam'
sam2=args.samPath + '/' + args.G1 + '_' + args.fq + '_upd_feature_uniq.sam'

sarray=[]
sarray.append(sam1)
sarray.append(sam2)

print(sarray)
### FQ files, num is 1 for SE, represents fq per 2 sam files

### Before Sam compare
with open(args.out, 'w') as outfile:

    if args.alnType=='SE':
        num=1
	if len(sarray) != 2*num:
            outfile.write(fq_file + ',' + 'Do NOT have 2 SAM files - rerun updated alignments!')
        else:
            #Continue with sam compare
            outfile.write(fq_file + ',' + 'Have 2 SAM files - good!')
    elif args.alnType=='PE':
        num=2
	if len(sarray)*2 != 2*num:
            outfile.write(fq_file + ',' + 'Do NOT have 2 SAM files - rerun updated alignments!')
        else:
            #Continue with sam compare
            outfile.write('Have 2 SAM files - good!')
    else:
	outfile.write('Select SE or PE!')
