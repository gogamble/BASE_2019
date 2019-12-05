#!/bin/sh
#SBATCH --mail-user=breccam@ufl.edu
#SBATCH --job-name=alnSim
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH --qos=mcintyre-b
#SBATCH -o /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/SLURM_LOGS/alnSim.%A-%a.out
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --mem=12gb
#SBATCH --array=84-273%10


#Set the number of runs that each SLURM task should do
PER_TASK=1

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo -e "This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM"

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ ))
do


    ## df is list of G2-G1-readLength combinations - 273 combinations
    # align simulated F1 reads to updated references
    # run sam compare script

    # Load modules
    module load seqkit/0.10.2
    module load bwa/0.7.7
    module load python/2.7.6
    module load bedtools/2.27.1
    module load samtools/1.9
    module load htslib

    ## Set / Create Directories and Variables
    PROJ=/ufrc/mcintyre/share/cegs2_MvsF_exp
    PY=$PROJ/scripts/ase_cegs_scripts/python_scripts
    SCRIPTS=$PROJ/scripts/ase_cegs_scripts
    REF=$PROJ/updated_genomes


    ## Assign paths to other input/output directories
    SIM=$PROJ/simulation
        mkdir -p $SIM
    SAMC=$SIM/ase_counts_updated_feature
        mkdir -p $SAMC
    CHKALN=$SIM/check_aln
        mkdir -p $SAMC
    CHKSC=$SIM/check_samcomp
        mkdir -p $CHKSC


    ## Get genotype info and fq files from design file, POSITIONAL
    DESIGN_FILE=$PROJ/design_files/df_g1_g2_RL_noHeader.csv
    DESIGN=$(sed -n "${RUN}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    G1=${ARRAY[0]}
    G2=${ARRAY[1]}
    READLEN=${ARRAY[2]}

    ## create temp directory
    ROZ=$PROJ/roz_sim_${RUN}
	mkdir -p $ROZ
    OUTALN=$ROZ/aln_upd_feature
        mkdir -p $OUTALN
    PSAM=$ROZ/bwa_parse
        mkdir -p $PSAM

    ### Cat together G1 and G2 simulated reads --> "F1"    (can take time)
    echo "
        catting G1 and G2 simulated reads
        "
    cat $SIM/simulated_reads/${G1}_${READLEN}bp_sim.fq $SIM/simulated_reads/${G2}_${READLEN}bp_sim.fq > $ROZ/sim_F1_${G2}_${G1}_${READLEN}bp.fq

    ## split into 20 files 
    echo "
	splitting into 20 separate files
        "
    mkdir -p $ROZ/split
    seqkit split2 $ROZ/sim_F1_${G2}_${G1}_${READLEN}bp.fq -p 20 -O $ROZ/split -f 

    echo "Remove catted file
    "
    rm $ROZ/sim_F1_${G2}_${G1}_${READLEN}bp.fq

    ###### (3) Align simulated split reads to updated genomes - first to G2 ref then to G1 reference - and Parse sam file
    for i in {001..020}
    do 
    
      READ=$ROZ/split/sim_F1_${G2}_${G1}_${READLEN}bp.part_${i}.fq
      FQ=sim_F1_${G2}_${G1}_${READLEN}bp.part_${i}
      FQFILE=${FQ}.fq
      echo " aligning ${READ}
           "

      for FOO in G1 G2
      do
          if [[ $FOO == 'G2' ]]
          then
              BREF=$REF/${G2}_snp_upd_genome_BWA
              echo $NUMREAD | awk -v fq=${FQFILE} -v gq=pre_aln_read_count '{print "filename" "," gq "\n" fq "," $0}' > $CHKALN/pre_aln_reads_${G2}_${FQ}.csv

              bwa mem -t 8 -M $BREF $READ > $ROZ/${G2}_${FQ}_upd_feature.sam

              module unload python/2.7.6
              module load python3
              $PY/BWASplitSAM_08amp.py -s $ROZ/${G2}_${FQ}_upd_feature.sam --outdir $PSAM -fq1 $READ

              ## cat together mapped and opposite
              cat $PSAM/${G2}_${FQ}_upd_feature_mapped.sam $PSAM/${G2}_${FQ}_upd_feature_opposite.sam > $PSAM/${G2}_${FQ}_upd_feature_uniq.sam


          elif [[ $FOO == 'G1' ]]
          then
              BREF=$REF/${G1}_snp_upd_genome_BWA
              echo $NUMREAD | awk -v fq=${FQ} -v gq=pre_aln_read_count '{print "filename" "," gq "\n" fq "," $0}' > $CHKALN/pre_aln_reads_${G1}_${FQ}.csv

              bwa mem -t 8 -M $BREF $READ > $ROZ/${G1}_${FQ}_upd_feature.sam
            
              module unload python/2.7.6
              module load python3
              $PY/BWASplitSAM_08amp.py -s $ROZ/${G1}_${FQ}_upd_feature.sam --outdir $PSAM -fq1 $READ

           ## cat together mapped and opposite
              cat $PSAM/${G1}_${FQ}_upd_feature_mapped.sam $PSAM/${G1}_${FQ}_upd_feature_opposite.sam > $PSAM/${G1}_${FQ}_upd_feature_uniq.sam

          fi
      done
      
      echo "Running python checks for ${FQ}
           "
      ### for every FQ file run, should have 2 sam files
      ## python script to check -
          module unload python3
          module load python/2.7.6

          python $SCRIPTS/check_sam_present_03brm.py \
            -fq $FQFILE \
            -alnType SE \
            -samPath $PSAM \
            -G1 $G1 \
            -G2 $G2 \
            -o $CHKALN/check_sam_present_${FQ}.txt

      ## run python script to count reads in (start reads) and reads out (end reads)
          python $SCRIPTS/check_for_lost_reads_04brm.py \
            -a1 $PSAM/${G1}_${FQ}_upd_feature_summary.csv \
            -a2 $PSAM/${G2}_${FQ}_upd_feature_summary.csv \
            -numread1 $CHKALN/pre_aln_reads_${G1}_${FQ}.csv \
            -numread2 $CHKALN/pre_aln_reads_${G2}_${FQ}.csv \
            -fq $FQFILE \
            -o $CHKALN/check_for_lost_reads_${FQ}.csv

      ######  Create modified BED file to use in SAM compare - this bed has features 1st with start and end positions in the chromosome
      ## use awk to reorder columms in bed file
      echo "Creating modified BED file for SAM compare
           "
      ## this bed file for SNPs only
      BEDFILE=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_fusions_coverage.bed

      SBED=$ROZ/snp_feature_first.bed

      awk -v OFS='\t' '{print $4,$2,$3,$1}' $BEDFILE > $SBED

      ###### Bedtools Intersect:   Here we will call the shell script to reformat the sam file so that the have feature names instead of CHR names
        ## In parsed SAM, 
        ## Sz158_sim_Sz158_w501_150bp_upd_feature_uniq.sam
        ## w501_sim_Sz158_w501_150bp_upd_feature_uniq.sam

      for SAMFILE in $PSAM/*_${FQ}_upd_feature_uniq.sam
      do
        MYSAMFILE2=$(basename $SAMFILE)

        AWKTMP=$PSAM/${MYSAMFILE2/_feature_uniq.sam/_feature_uniq_AWK.txt}
        NEWSAM=$PSAM/${MYSAMFILE2/_feature_uniq.sam/_feature_uniq_REAL.sam}

        #Create a bed file to write the  starting position of every read
        awk -v OFS='\t' '{print $3,$4,$4}' $SAMFILE > $AWKTMP
        BED3=$PSAM/${MYSAMFILE2/_feature_uniq.sam/_feature_uniq_int.bed}

        #Run bedtools intersect with -loj between the reads and the features.
        #We will have one result for each region
        bedtools intersect -a $AWKTMP -b $BEDFILE -loj > $BED3

        #With awk substitute column 3 of sam file with column 7 (Feature name) of bed file.
        awk -v OFS='\t' 'FNR==NR{a[NR]=$7;next}{$3=a[FNR]}1' $BED3 $SAMFILE | awk '$3!="."' > $NEWSAM

        echo initial sam file $SAMFILE
        echo awk outfile $AWKTMP
        echo bed intersect outfile $BED3
        echo new sam file "$NEWSAM"

      done

      ###### Run Sam Compare
      # grab sam files

      SAM1=$PSAM/${G1}_${FQ}_upd_feature_uniq_REAL.sam
      SAM2=$PSAM/${G2}_${FQ}_upd_feature_uniq_REAL.sam

      READ1=$FQFILE

      echo -e "READ1: '${READ1}"
      echo -e "SAM1: '${SAM1}'"
      echo -e "SAM2: '${SAM2}'"
      echo -e "BED: '${SBED}'
           "

      echo "Starting sam compare for $G2 and $READ
           "
      python $PY/sam_compare_w_feature.py \
        -n \
        -l ${READLEN} \
        -f $SBED \
        -q $READ1 \
        -A $SAM1 \
        -B $SAM2 \
        -c $SAMC/ase_counts_${FQ}.csv \
        -t $SAMC/ase_totals_${FQ}.txt \
        --log $CHKSC/ase_log_${FQ}.log
    
    done

   ## delete temp dir
   rm -r $ROZ

done
