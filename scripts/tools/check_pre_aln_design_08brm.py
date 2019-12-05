import argparse, csv
import pandas as pd
import os
import sys

def getOptions():
    parser = argparse.ArgumentParser(description='Run checks on design file as prep for ase_counts_sum')
    parser.add_argument('-design','--design',dest='design', action='store', required=True, help='Design file containing fq file names and sample ids [Required]')
    parser.add_argument('-id','--id',dest='id', action='store', required=True, help='Name of the column containing sampleIDs [Required]')
    parser.add_argument('-fnames','--fnames',dest='fnames', action='store', required=True, help='Name of the column containing fastq file names [Required]')
    parser.add_argument('-g1','--g1',dest='g1', action='store', required=True, help='Name of the column containing G1 names [Required]')
    parser.add_argument('-g2','--g2',dest='g2', action='store', required=True, help='Name of the column containing G2 names [Required]')
    parser.add_argument('-rep','--rep',dest='rep', action='store', required=True, help='fastq extension (ex: .fq, .fastq) [Required]')
    parser.add_argument('-ex','--ex',dest='ex', action='store', required=True, help='Name of the column containing rep values [Required]')
    parser.add_argument('-dir','--dir',dest='dir', action='store', required=True, help='Give the entire path to the directory where design and fastq files located [Required]')
    parser.add_argument('-readlen','--readlen',dest='readlen', action='store', required=True, help='Name of column containing readLength values [Required]')
    args = parser.parse_args()
    return(args)

#check that fastq files are listed uniquely and exist in directory; returns duplicate rows if they exist
def fastq_check(design, file_col, path, ext):
    df = pd.read_csv(os.path.join(path,design))
    if len(df[file_col].unique().tolist()) < len(df[file_col].tolist()): ##creates lists out of elements in column and compares to a list of unique elements
        dups = df[df.duplicated(file_col, keep=False) == True] ## creates list of any duplicated rows
        print("Error: duplicate rows:")
        print(dups)
      	sys.exit()

    count2 = 0
    for i in range(len(df)): ## row by row iteration
        filename = df.iloc[i].loc[file_col] ## get filename from column containing filenames
        file = filename + ext ## append user-provided extension to filename
        gz = filename + '.gz' ## create gz variable to check that files do not exist as gz
        if os.path.isfile(os.path.join(path, file)):
            count2 += 1
            continue 
        else:
            if os.path.isfile(os.path.join(path, gz)):
                print("Error: file " + gz + " is zipped. Please unzip your file") ##checks if file is zipped (gz) and notifies user
                sys.exit()
            else:
                print("File " + file + " does not exist in directory " + path)
                sys.exit()
		
    print("rows checked2:", count2)

#check that columns G1, G2, sampleName, comparate_1, comparate_2, rep, readLength, and fqName exist in design file (and in that order!)
#Cannot check if user mislabeled column
def columns_check(design, g1_col, g2_col, sampleID_col, file_col, rep_col, readlen_col, path):
    df = pd.read_csv(os.path.join(path,design))
    if g1_col not in df:
        print("Error: column " + g1_col + " does not exist in design file")
        sys.exit()
    if g2_col not in df:
        print("Error: column " + g2_col + " does not exist in design file")
        sys.exit()
    if sampleID_col not in df: 
        print("Error: column " + sampleID_col + " does not exist in design file")
        sys.exit()
    if file_col not in df:
        print("Error: column " + file_col + " does not exist in design file")
        sys.exit()
    if rep_col not in df:
        print("Error: column " + rep_col + " does not exist in design file")
        sys.exit()
    if readlen_col not in df:
        print("Error: column " + readlen_col + " does not exist in design file")
        sys.exit()
    columns = ['G1','G2','sampleID','fqName','fqExtension','techRep','readLength']
    headers = list(df)
    if headers != columns:
        print("ERROR: column headers in file " + design + " do not align with order requirements, please check.")
        sys.exit()

def main(args):
    columns_check(args.design, args.g1, args.g2, args.id, args.fnames, args.rep, args.readlen, args.dir)
    fastq_check(args.design, args.fnames, args.dir, args.ex)

##run main
args = getOptions()
main(args)
