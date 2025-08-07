#!/bin/bash
#SBATCH --job-name="salix_plastome_phylogeny"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --array=52,53,62,66,68
#SBATCH --mem=1G #this is per job

function index {

REFDIR=$1

cd $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}
cat $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/assemble_out/*fasta > $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/tmp.fasta
bowtie2-build -f $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/tmp.fasta ${prefix}_${SLURM_ARRAY_TASK_ID}
rm $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/tmp.fasta
}

function map {

INDIR=$1
REFDIR=$2

cd $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}
#bowtie2-align-s --wrapper basic-0 -x ${prefix}_${SLURM_ARRAY_TASK_ID} -p 20 --threads 16 --phred33 --very-sensitive --quiet --time --rg-id ${prefix} --rg SM:${prefix} --rg PL:'ILLUMINA' -S $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/${prefix}_${SLURM_ARRAY_TASK_ID}.bam -1 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz -2 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz
samtools view -Sbh -F 4 $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/${prefix}_${SLURM_ARRAY_TASK_ID}.bam|bamtools filter -tag XM:'<=1'|samtools sort -o $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
samtools index $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
samtools coverage $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam > $REFDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/${prefix}_${SLURM_ARRAY_TASK_ID}_depth.txt
}

function main {

#UNIVERSAL VARIABLES
env_name=easy353
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_combined_salix #the last time i just ran it with salix_47 by mistake
prefix=salix
READS=$HOME/projects/rbge/zedchen/barcoding/salix_combined/results/00_reads
#SET UP DIRECTORIES
easy353path=$HOME/projects/rbge/$USER/env/$env_name/bin/Easy353
PARENT=$HOME/scratch/$WORKDIR/results
ref353=$PARENT/353_ref_Salix
#subdirectories
RESULT1=$SCRATCH/$WORKDIR/results/01_easy353_compile
RESULT2=$SCRATCH/$WORKDIR/results/02_easy353_alignment
RESULT3=$SCRATCH/$WORKDIR/results/03_easy353_singlegene_phy
RESULT4=$SCRATCH/$WORKDIR/results/04_easy353_astral
RESULT5=$SCRATCH/$WORKDIR/results/05_mapping_depth_check
#IN ARRAY
#index $RESULT5 $RESULT5 
map $READS $RESULT5 
}


#execution
#SLURM_ARRAY_TASK_ID=52
main
