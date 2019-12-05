#!/bin/sh
#SBATCH --mail-user=breccam@ufl.edu
#SBATCH --job-name=sim
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/SLURM_LOGS/sim.%A-%a.out
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=500M
#SBATCH --array=2-11,50-93

## simulate G1 and G2 reads from updated genomes!
## doing 50, 60 and 150 bp read lengths since this is what actual data consists of


## Load modules
module load python/2.7.6

## Set / Create Directories and Variables
PROJ=/ufrc/mcintyre/share/cegs2_MvsF_exp
PY=$PROJ/scripts/python_scripts
REF=$PROJ/updated_genomes
SCRIPTS=$PROJ/scripts/ase_cegs_scripts
SIM=$PROJ/simulation
    mkdir -p $SIM
OUT=$SIM/simulated_reads
    mkdir -p $OUT 

#### simulate reads 

## Get Genomes for simulating reads from design file
DESIGN_FILE=$PROJ/design_files/df_ase_cegs_genome_list.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

GENOME=${ARRAY[0]}
READLEN="50 60 150"
for i in $READLEN 
do

    $SCRIPTS/simulate_reads.py \
        --input $REF/${GENOME}_snp_upd_genome.fasta \
        --out $OUT/${GENOME}_${i}bp_sim.fq \
        --num ${i} \
        --prefix ${GENOME}${i}bp
done
