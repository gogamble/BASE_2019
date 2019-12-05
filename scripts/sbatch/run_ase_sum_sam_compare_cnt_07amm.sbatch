#!/bin/sh
#SBATCH --mail-user=ammorse@fl.edu
#SBATCH --job-name=sumSam
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/SLURM_LOGS/sumSam.%A.out
#SBATCH -t 6:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=6G


# design file MUST have column called 'sampleID' and a column called 'fqName'
# for each sampleID, sum the ase_counts for all 'tech reps' (sum different lanes, read lengths etc for same sample) 
# do NOT sum biological replicates
# A separate step because all count tables must be made before starting

# Load modules
module load python3
module load htslib

## Set / Create Directories and Variables
PROJ=/ufrc/mcintyre/share/cegs2_MvsF_exp
SAMC=$PROJ/simulation/ase_counts_updated_feature
SCRIPTS=$PROJ/scripts/ase_cegs_scripts
BED=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_fusions_coverage.bed

## Create output directory for summed sam count tables
SSUM=$SAMC/ase_counts_summed
    if [ ! -e $SSUM ]; then mkdir -p $SSUM; fi

## Indicate whether count tables to be summed are created from simulated data or original data (True=Simulated or False=Original Data)
SIM=True

## Set design file paths
DF=$PROJ/design_files
#DESIGN=$DF/cegs_MF_df_from_user_sorted.csv
DESIGN=$DF/df_sum_sim_count_tables.csv
## Set subset start and end positions (design files have headers, start at 2)

STARTRUN=2
RUNSIZE=5441

echo "
    summing ase counts for all genotypes
    "
## running all, on data
python3 $SCRIPTS/combine_cnt_tables_12jeb.py \
    -design $DESIGN \
    -sim $SIM \
    -bed $BED \
    -path $SAMC \
    -designdir $DF \
    -out $SSUM
