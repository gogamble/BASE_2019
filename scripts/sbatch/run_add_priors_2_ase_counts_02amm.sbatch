#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=sumSamSim
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/SLURM_LOGS/sumSamSim_%A-%a.out
#SBATCH -t 06:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=6gb
#SBATCH --array=2,17,18,46,79

### DF has header!! 

### (1) calculate priors from simulation

# design file MUST have column called 'sampleID' and a column called 'fqName' and a column called 'readLength'

# Load modules
module purge
#module load python3
module load cuda/9.2.88   # needed to load python3 version correctly
module load intel/2018.1.163 # needed to load python3 version correctly
module load python3/3.6.2019.4.088
module load htslib


## Set / Create Directories and Variables
PROJ=/ufrc/mcintyre/share/cegs2_MvsF_exp
SCRIPTS=$PROJ/scripts/ase_cegs_scripts
SIM=$PROJ/simulation
SIMC=$SIM/ase_counts_updated_feature
SSUM=$SIMC/ase_counts_summed

PRIOR=$PROJ/simulation/priors_02amm
    mkdir -p $PRIOR

## Set design file directory
DESIGN_FILE=$SSUM/df_sim_samcomp_summed.csv
DESIGN=$(cat $DESIGN_FILE | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

G1=${ARRAY[0]}
G2=${ARRAY[1]}

##### (1) calculate priors
python3 $SCRIPTS/calculate_priors_loop.py \
    --output $PRIOR \
    --g1 $G1 \
    --g2 $G2 \
    --sam-compare-dir $SSUM  \
    --debug
