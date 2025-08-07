#!/bin/bash
#SBATCH --job-name="salix_plastome_phylogeny"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --array=1-28 
#SBATCH --mem=32G #change it to the amount of memory you need

function make_work_dir { 

#mkdir $PARENT

#Create folders under result folder for each program
#QC results


mkdir $RENAMED
mkdir $RESULT1
mkdir $RESULT2
mkdir $RESULT3

}


#Whole plastome assembly: array
function plastome_assembly {

INDIR=$1
OUTDIR=$2
ref=$3

get_organelle_from_reads.py -s $ref -1 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz -2 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
                            -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID} -R 15 -k 21,45,65,85,105 -F embplant_pt

}


function main {

#UNIVERSAL VARIABLES
env_name=salix
USER=zedchen 
WORKDIR=Barcoding_km
prefix=salix
#SET UP DIRECTORIES
adapters=/home/zchen/apps/conda/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa #need to specify the full path to the adapter files
DATA=$HOME/scratch/$WORKDIR/data/willow_seq
ref=$HOME/scratch/$WORKDIR/reference/Salix_blinii_plastome.fasta
PARENT=$HOME/scratch/$WORKDIR/results
#subdirectories
RENAMED=$HOME/scratch/$WORKDIR/results/renamed
RESULT1=$HOME/scratch/$WORKDIR/results/qc1_trimmomatic #this step is omitted eventually
RESULT2=$HOME/scratch/$WORKDIR/results/qc2_fastp
RESULT3=$HOME/scratch/$WORKDIR/results/phy1_plastome
RESULT4=$HOME/scratch/$WORKDIR/results/phy2_alignment_tree

#make_work_dir 
#rename the files for submitting in slurm array: run only once, the copying process takes a while

#IN ARRAY
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >> plastome_debug.log
plastome_assembly $RESULT2 $RESULT3 $ref #you have to write the inputs

}


#execution
main