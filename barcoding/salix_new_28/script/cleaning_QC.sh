#!/bin/bash
#SBATCH --job-name="salix_barcoding_QC"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --array=5-28 
#SBATCH --mem=32G #for each task. not in total

function make_work_dir { 

#mkdir $PARENT

#Create folders under result folder for each program
#QC results


mkdir $RENAMED
mkdir $RESULT1
mkdir $RESULT2

}

function trim_adapter {

INDIR=$1
OUTDIR=$2

trimmomatic PE ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
                ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
                ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

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
WORKDIR=Barcoding_km/barcoding_salix_28
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
#write a csv file for matching the names back later
#./rename_files.py -i $DATA -o $RENAMED --prefix $prefix

#IN ARRAY
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >> debug.log
#trim_adapter $RENAMED $RESULT1
fastp_qc $RENAMED $RESULT2 #you have to write the inputs
#
}

#UNIVERSAL VARIABLES


#execution
main