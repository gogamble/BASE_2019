#!/bin/sh
#SBATCH --mail-user=breccam@ufl.edu
#SBATCH --job-name=mergeComp
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/SLURM_LOGS/mergeComp_%A-%a.out
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=2gb
#SBATCH --array=2

# Load modules
module load python3
module load htslib



## Set / Create Directories and Variables
PROJ=/ufrc/mcintyre/share/cegs2_MvsF_exp
SCRIPTS=$PROJ/scripts/ase_cegs_scripts
SIM=$PROJ/simulation

BAYESIN=$SIM/bayesian_in
COMPIN=$BAYESIN/comparison_input
    if [ ! -e $COMPIN ]; then mkdir -p $COMPIN; fi

DF=$PROJ/ase_counts_updated_feature/ase_counts_summed_02amm
DESIGN_FILE=$DF/df_ase_for_bayesian_modified.csv

## Set subset start and end positions (design files have headers, start at 2)

STARTRUN=2
RUNSIZE=5

#### (1) Read user-inputted comparisons and merge corresponding comparate files to be analyzed by the Bayesian Model


python3 $SCRIPTS/merge_comparates_for_bayesian_03amm.py \
    --output $COMPIN \
    --comp $BAYESIN \
    --design $DESIGN_FILE \
    --begin $STARTRUN \
    --end $RUNSIZE 
