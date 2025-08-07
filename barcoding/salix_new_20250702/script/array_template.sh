#!/bin/bash
#SBATCH --job-name="barcoding_QC"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --array=1-28 
#SBATCH --mem=128G #change it to the amount of memory you need

function make_work_dir { 

#mkdir $PARENT

#Create folders under result folder for each program
#QC results


mkdir $RENAMED
mkdir $RESULT1
mkdir $RESULT2

}



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

function main {

env_name=salix
USER=zedchen 
WORKDIR=Barcoding_km
prefix=salix
#SET UP DIRECTORIES
adapters=/home/zchen/apps/conda/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa #need to specify the full path to the adapter files
DATA=$HOME/scratch/$WORKDIR/data/willow_seq
PARENT=$HOME/scratch/$WORKDIR/results
#subdirectories
RENAMED=$HOME/scratch/$WORKDIR/results/renamed
RESULT1=$HOME/scratch/$WORKDIR/results/qc1_trimmomatic #this step is omitted eventually
RESULT2=$HOME/scratch/$WORKDIR/results/qc2_fastp

#make_work_dir 
#rename the files for submitting in slurm array: run only once, the copying process takes a while

#IN ARRAY
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >> debug.log
fastp_qc $RENAMED $RESULT2 #you have to write the inputs
#
}

#UNIVERSAL VARIABLES


#execution
main