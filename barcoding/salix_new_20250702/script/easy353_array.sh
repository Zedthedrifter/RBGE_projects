#!/bin/bash
#SBATCH --job-name="salix47easy353"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --array=1-13
#SBATCH --partition=short
#SBATCH --mem=32G #change it to the amount of memory you need

function fastp_qc {

#SKIP TRIM ADAPTOR
INDIR=$1
OUTDIR=$2

fastp -i ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz -I ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
       -o ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz -O ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
       --unpaired1 ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_unpaired1.fq.gz --unpaired2 ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_unpaired2.fq.gz \
       -j ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.json \
       -h ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.html \
       --dont_overwrite --detect_adapter_for_pe -c > ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.fastp.log


}

function easy353_assembly {

INDIR=$1
OUTDIR=$2
ref353=$3

$easy353path/easy353.py -1 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz -2 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
                        -r $ref353/353gene -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID} \
                        -fk 31 -ak 41 -ft 1 -at 4 -kmer_limit 8 


}

function main {

#UNIVERSAL VARIABLES
env_name=easy353
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_salix_13
prefix=salix
#SET UP DIRECTORIES
easy353path=$HOME/projects/rbge/$USER/env/$env_name/bin/Easy353
PARENT=$HOME/scratch/$WORKDIR/results
ref353=$PARENT/353_ref_Salix
#subdirectories
RENAMED=$HOME/scratch/$WORKDIR/results/renamed
RESULT1=$HOME/scratch/$WORKDIR/results/qc1_trimmomatic #this step is omitted eventually
RESULT2=$HOME/scratch/$WORKDIR/results/qc2_fastp
RESULT3=$HOME/scratch/$WORKDIR/results/phy1_plastome
RESULT4=$HOME/scratch/$WORKDIR/results/phy2_alignment_tree
RESULT5=$HOME/scratch/$WORKDIR/results/phy3_easy353assembly


#fastp_qc $RENAMED $RESULT2 #
easy353_assembly $RESULT2 $RESULT5 $ref353 #use fastp cleaned reads

}


#execution
main