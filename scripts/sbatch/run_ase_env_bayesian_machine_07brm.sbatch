#!/bin/sh
#SBATCH --mail-user=breccam@ufl.edu
#SBATCH --job-name=env_bayesian
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH --qos=mcintyre-b
#SBATCH -o /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/SLURM_LOGS/env_bayesian.%A-%a.out
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --array=5 

# Design file MUST have a column for "G1" "G2" comparisons, and "compID"
# G1 G2 cannot start with a number or special character
# MUST be alphaneumeric, and "L" "T" "M" "V" "line" "tester" cannot be used due to original code constraints

# Load modules
module purge
#module load R/3.4.3
#module load python3
module load cuda/9.2.88   # needed to load python3 version correctly
module load intel/2018.1.163 # needed to load python3 version correctly
module load python3/3.6.2019.4.088
module load htslib
module load R/3.6

## Set / Create Directories and Variables
PROJ=/ufrc/mcintyre/share/cegs2_MvsF_exp
SCRIPTS=$PROJ/scripts/ase_cegs_scripts

## Create temp output directory for initial datafile with modified headers before entering Bayesian
ROZ=$PROJ/roz_test
    mkdir -p $ROZ

## Create output directory for AI data
ENV=$PROJ/bayesian_out/launcher_model
    if [ ! -e $ENV ]; then mkdir -p $ENV; fi

## Set Bayesian subprocess path (equivalent to model_launcher.r), MUST BE IN QUOTATIONS
SUBPATH='/ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/ase_cegs_scripts/R_scripts/model_launcher.r'
#SUBPATH='/ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/ase_cegs_scripts/simul_NBmodel_WITHBOTH8_gammaalphadeltainboth_taugamma1_withSTAN4_02brm_ncond.R'
#SUBPATH='/ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/ase_cegs_scripts/simul_NBmodel_WITHBOTH8_gammaalphadeltainboth_taugamma1_withSTAN4_03brm_ncond.R'
#SUBPATH='/ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/ase_cegs_scripts/2019_stan_V04/simul_NBmodel_WITHBOTH8_gammaalphadeltainboth_taugamma1_withSTAN4_02brm_3cond.R'

#SUBPATH='/ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/ase_cegs_scripts/2019_stan_V04/simul_NBmodel_WITHBOTH8_gammaalphadeltainboth_taugamma1_withSTAN4_02brm_3cond.R'
## Set design file with G1, G2, c1, c2, and input filename columns
DESIGN_FILE=$PROJ/ase_counts_updated_feature/ase_counts_summed_02amm/df_ase_for_bayesian_modified.csv

#DESIGN_FILE=$PROJ/ase_counts_updated_feature/ase_counts_summed_02amm/df_ase_for_bayesian_rest.csv
#DESIGN_FILE=$PROJ/ase_counts_updated_feature/ase_counts_summed_02amm/df_ase_for_bayesian.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

ARRAY=$SLURM_ARRAY_TASK_ID
echo $ARRAY

DATAIN=$PROJ/simulation/bayesian_in/comparison_input


###### (1) Check datafile columns, g1 g2  c1 c2 and filename inputs
#python3 $SCRIPTS/check_env_bayesian_design_05brm.py \
#    -design $DESIGN_FILE \
#    -dir $PROJ/design_files/

###### (2) Run python script calling environmental bayesian model
python3 $SCRIPTS/NBmodel_16brm.py \
    -datafile $DATAIN/ \
    -datafile2 $ROZ/ \
    -design $DESIGN_FILE \
    -routput $ENV/ \
    -subpath $SUBPATH \
    -array $ARRAY \
    -o $ENV/ 
