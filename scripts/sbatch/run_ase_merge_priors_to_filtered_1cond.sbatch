#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=priorfilt
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/SLURM_LOGS/priorfilt_%A-%a.out
#SBATCH -t 06:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=6gb
#SBATCH --array=2-471

# Load modules
module load python3
module load htslib



## Set / Create Directories and Variables
PROJ=/ufrc/mcintyre/share/cegs2_MvsF_exp
SCRIPTS=$PROJ/scripts/ase_cegs_scripts
#ROZ=$PROJ/roz_sim
#SAMC=$PROJ/simulation/ase_counts_updated_feature
SIM=$PROJ/simulation

SSUM=$PROJ/ase_counts_summarized

DF=$PROJ/ase_counts_updated_feature/ase_counts_summed_02amm
DESIGN_FILE=$DF/df_ase_samcomp_summed_1cond.csv

DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

G1=${ARRAY[0]}
G2=${ARRAY[1]}
COMP=${ARRAY[2]}
echo $COMP

## TO DO:
## Make automated design file

#### (1) merge summarized ase counts to priors
BAYESIN=$SIM/bayesian_in
    mkdir -p $BAYESIN
PRIOR=$SIM/priors

## Remember to switch output back to bayesin

python3 $SCRIPTS/merge_sim_filtered_1cond_02brm.py \
    --output $BAYESIN \
    --countsf $SSUM \
    --comp $COMP \
    --sf $PRIOR \
    --g1 $G1 \
    --g2 $G2
