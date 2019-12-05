#!/bin/sh
#SBATCH --mail-user=breccam@ufl.edu
#SBATCH --job-name=filt
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/SLURM_LOGS/filt.%A.out
#SBATCH -t 11:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=6gb


#### filter ase counts tables where techReps have been summed
### ase_compare column must have "rep" or numerical indication removed after last underscore
## DF:
	## G1,G2,sampleID,ase_compare
	## w501,Sz11,Sz11_M_etoh_time10_rep1,Sz11_F_noetoh_time10

# Load modules
module load python3
module load htslib
module load R

## Set directories / variables
PROJ=/ufrc/mcintyre/share/cegs2_MvsF_exp
SCRIPTS=$PROJ/scripts/ase_cegs_scripts
SSUM=$PROJ/ase_counts_updated_feature/ase_counts_summed_02amm

DESIGN=$SSUM/df_ase_samcomp_summed_1cond.csv
APN=1

## Create output directory for filtered ase counts
FILT=$PROJ/ase_counts_summarized
    mkdir -p $FILT

#TO DO:
## NEED TO ADD SUBSET OPTION

python3 $SCRIPTS/summarize_sam_compare_cnts_table_1cond_03brm.py \
    --output $FILT \
    --design $DESIGN \
    --compCol Comparate2 \
    --lineCol G2 \
    --sam-compare-dir $SSUM \
    --apn $APN \
