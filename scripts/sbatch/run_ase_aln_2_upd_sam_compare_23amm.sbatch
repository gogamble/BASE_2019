#!/bin/sh
#SBATCH --mail-user=breccam@ufl.edu
#SBATCH --job-name=aln_upd
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/SLURM_LOGS/aln_upd.%A-%a.out
#SBATCH -t 6:00:00
#SBATCH --ntasks-per-node=8
#SBATCH --mem=2gb
###SBATCH --array=1-423
#SBATCH --array=1


# aligns to updated references using bash loop across design file
# runs sam compare script

# Load modules
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
ORIG=$PROJ/cutadapt_bbmerge_fastqSplitDups
BEDFILE=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_fusions_coverage.bed


## Create output directory for alignments to updated_feature
OUTALN=$PROJ/aln_upd_feature
    if [ ! -e $OUTALN ]; then mkdir -p $OUTALN; fi

PSAM=$OUTALN/bwa_parse
    if [ ! -e $PSAM ]; then mkdir -p $PSAM; fi

## Create output directory for ase counts from sam compare
SAMC=$PROJ/ase_counts_updated_feature
    if [ ! -e $SAMC ]; then mkdir -p $SAMC; fi

## Create output diretory for alignment summaries
ALNSUM=$OUTALN/align_summaries
    mkdir -p $ALNSUM

## Create output directory for upd aln checks
CHKALN=$OUTALN/check_aln
    if [ ! -e $CHKALN ]; then mkdir -p $CHKALN; fi

## Create output directory for Sam Compare check
CHKSC=$SAMC/check_samcomp
    if [ ! -e $CHKSC ]; then mkdir -p $CHKSC; fi


## Using looping of the array task runs
## df (cegs_MF_df_from_user_03brm_noHeader.csv) has 4230 rows
## each row takes 3 minutes to run so want 5 rows in the loop (3 min *10rows  = 30 min)
        ## split into X array tasks:  4230/X = 10 rows per array task, X=423 (ran out of time with per_task=10)

#Set the number of runs that each SLURM task should do
PER_TASK=10

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo -e "
        This is task arrayID $SLURM_ARRAY_TASK_ID, which will do rows $START_NUM to $END_NUM in the DF
         "

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ ))
do
    echo "
        we are on row: $RUN
        "
    DF=$PROJ/design_files
    DESIGN_FILE=$DF/cegs_MF_df_from_user_03brm_noHeader.csv
    DF_HEADER=$DF/cegs_MF_df_from_user_03brm.csv

    ## Check design file that it meets requirements --> THIS WONT RUN ON DESIGN FILE WITHOUT HEADER
        python $SCRIPTS/check_pre_aln_design_08brm.py \
            -design ${DF_HEADER} \
            -id sampleID \
            -fnames fqName \
            -g1 G1 \
            -g2 G2 \
            -rep techRep \
            -dir $ORIG \
            -ex .fq \
            -readlen readLength

    ## Get genotype info and fq files from design file, POSITIONAL
    DESIGN=$(sed -n "${RUN}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    G1=${ARRAY[0]}
    G2=${ARRAY[1]}
    READLEN=${ARRAY[6]}
    FQ=${ARRAY[3]}

    ## create temp directory 
    ROZ=$PROJ/roz_test_${RUN}_${FQ}
        mkdir -p $ROZ

    SBED=$ROZ/snp_feature_first.bed

    ######  Create modified BED file to use in SAM compare - this bed has features 1st with start and end positions in t$
    ## use awk to reorder columms in bed file
    echo "
        Reformatting $BEDFILE"
        awk -v OFS='\t' '{print $4,$2,$3,$1}' $BEDFILE > $SBED


    ###### (1) Align Reads to Updated Genomes - first to G2 ref then to G1 reference - and Parse sam file
    READ=$ORIG/${FQ}.fq
    FQLINEFN=$(wc -l $READ)
    FQLINE=$(echo $FQLINEFN | cut -d" " -f1)
    NUMREAD=$(( FQLINE / 4 )) 
    FN=$(echo $FQLINEFN | cut -d" " -f2)

    for FOO in G1 G2
    do
        if [[ $FOO == 'G2' ]]
        then
            BREF=$REF/${G2}_snp_upd_genome_BWA

            echo $NUMREAD | awk -v fq=${FQ} -v gq=pre_aln_read_count '{print "filename" "," gq "\n" fq "," $0}' > $CHKALN/pre_aln_reads_${G2}_${FQ}.csv

            echo "
                 Aligning to output  $OUTALN/${G2}_${FQ}_upd_feature.sam"
            bwa mem -t 8 -M $BREF $READ > $OUTALN/${G2}_${FQ}_upd_feature.sam

            echo "
                 Start BWASplitSam  $OUTALN/${G2}_${FQ}_upd_feature.sam"
            module unload python/2.7.6
            module load python3
            $PY/BWASplitSAM_08amp.py -s $OUTALN/${G2}_${FQ}_upd_feature.sam --outdir $PSAM -fq1 $READ
                 
            ## move alignment summaries to new directory

            ## cat together mapped and opposite
            cat $PSAM/${G2}_${FQ}_upd_feature_mapped.sam $PSAM/${G2}_${FQ}_upd_feature_opposite.sam > $PSAM/${G2}_${FQ}_upd_feature_uniq.sam

        elif [[ $FOO == 'G1' ]]
        then
            BREF=$REF/${G1}_snp_upd_genome_BWA
            echo $NUMREAD | awk -v fq=${FQ} -v gq=pre_aln_read_count '{print "filename" "," gq "\n" fq "," $0}' > $CHKALN/pre_aln_reads_${G1}_${FQ}.csv
            echo "
                 Aligning to output  $OUTALN/${G1}_${FQ}_upd_feature.sam"
            bwa mem -t 8 -M $BREF $READ > $OUTALN/${G1}_${FQ}_upd_feature.sam

            echo "
                 Start BWASplitSam  $OUTALN/${G2}_${FQ}_upd_feature.sam"
            module unload python/2.7.6
            module load python3
            $PY/BWASplitSAM_08amp.py -s $OUTALN/${G1}_${FQ}_upd_feature.sam --outdir $PSAM -fq1 $READ

            ## cat together mapped and opposite
            cat $PSAM/${G1}_${FQ}_upd_feature_mapped.sam $PSAM/${G1}_${FQ}_upd_feature_opposite.sam > $PSAM/${G1}_${FQ}_upd_feature_uniq.sam
        fi
    done

    ## move alignment summaries to new directory
    mv $PSAM/*_summary.csv $ALNSUM/

    ### for every FQ file run, should have 2 sam files
    ## python script to check - 

#    module unload python3
#    module load python/2.7.6
#    python $SCRIPTS/check_sam_present_03brm.py \
#        -fq $FQ \
#        -alnType SE \
#        -samPath $PSAM \
#        -G1 $G1 \
#        -G2 $G2 \
#        -o $CHKALN/check_sam_present_${FQ}.txt

    ## run python script to count reads in (start reads) and reads out (end reads)

#    python $SCRIPTS/check_for_lost_reads_04brm.py \
#        -a1 $PSAM/${G1}_${FQ}_upd_feature_summary.csv \
#        -a2 $PSAM/${G2}_${FQ}_upd_feature_summary.csv \
#        -numread1 $CHKALN/pre_aln_reads_${G1}_${FQ}.csv \
#        -numread2 $CHKALN/pre_aln_reads_${G2}_${FQ}.csv \
#        -fq $FQ \
#        -o $CHKALN/check_for_lost_reads_${FQ}.csv


    ## Insert bedtools intersect script + any checks  (reads in aln sam output)
    ###### (2) Bedtools Intersect:   Here we will call the shell script to reformat the sam file so that the have feature names instead of CHR names
    ## In parsed SAM, 

#    for SAMFILE in $PSAM/*_${FQ}_upd_feature_uniq.sam
#    do
#        MYSAMFILE2=$(basename $SAMFILE)

#        AWKTMP=$PSAM/${MYSAMFILE2/_feature_uniq.sam/_feature_uniq_AWK.txt}
#        NEWSAM=$PSAM/${MYSAMFILE2/_feature_uniq.sam/_feature_uniq_REAL.sam}    

        #Create a bed file to write the  starting position of every read
#        awk -v OFS='\t' '{print $3,$4,$4}' $SAMFILE > $AWKTMP
#        BED3=$PSAM/${MYSAMFILE2/_feature_uniq.sam/_feature_uniq_int.bed}

        #Run bedtools intersect with -loj between the reads and the features. 
        #We will have one result for each region 
#        bedtools intersect -a $AWKTMP -b $BEDFILE -loj > $BED3
    
        #With awk substitute column 3 of sam file with column 7 (Feature name) of bed file.
#        awk -v OFS='\t' 'FNR==NR{a[NR]=$7;next}{$3=a[FNR]}1' $BED3 $SAMFILE | awk '$3!="."' > $NEWSAM
#        echo initial sam file $SAMFILE 
#        echo awk outfile $AWKTMP
#        echo bed intersect outfile $BED3
#        echo new sam file "$NEWSAM"

#    done 

       

    ## Grab sam files and bed files

#    SAM1=$PSAM/${G2}_${FQ}_upd_feature_uniq_REAL.sam
#    SAM2=$PSAM/${G1}_${FQ}_upd_feature_uniq_REAL.sam
#    BED1=$PSAM/${G2}_${FQ}_upd_feature_uniq_int.bed
#    BED2=$PSAM/${G1}_${FQ}_upd_feature_uniq_int.bed

#    awk 'NR==FNR{c[$3]++;next};c[$7] == 0' $SAM1 $BED1 > $CHKSC/check_sam_bed_${G2}_${FQ}.txt
#    awk 'NR==FNR{c[$3]++;next};c[$7] == 0' $SAM2 $BED2 > $CHKSC/check_sam_bed_${G1}_${FQ}.txt

    ###### (3) Run Sam Compare

#    READ1=${ORIG}/$FQ

#    echo -e "READ1: '${READ1}"
#    echo -e "SAM1: '${SAM1}'"
#    echo -e "SAM2: '${SAM2}'"
#    echo -e "BED: '${SBED}'"

#    echo " starting sam compare for $G2 and $FQ "
#    python $PY/sam_compare_w_feature.py \
#        -n \
#        -l ${READLEN} \
#        -f $SBED \
#        -q $READ1 \
#        -A $SAM1 \
#        -B $SAM2 \
#        -c $SAMC/ase_counts_${FQ}.csv \
#        -t $SAMC/ase_totals_${FQ}.txt \
#        --log $CHKSC/ase_log_${FQ}.log

#    echo " run sam compare check for $G2 and $FQ "
    # Check to make sure counts in csv summary file is within range of minimum unique reads from respective sam files and 
    # the summation of the unique reads of both sam files
#    python $SCRIPTS/check_samcomp_for_lost_reads_02brm.py \
#       -b1 $PSAM/${G1}_${FQ}_upd_feature_summary.csv \
#       -b2 $PSAM/${G2}_${FQ}_upd_feature_summary.csv \
#       -G1 $G1 \
#       -G2 $G2 \
#       -s $SAMC/ase_totals_${FQ}.txt \
#       -fq $FQ \
#       -o $CHKSC/check_samcomp_${FQ}_aln_2_upd.txt

    rm -r $ROZ
    
done
