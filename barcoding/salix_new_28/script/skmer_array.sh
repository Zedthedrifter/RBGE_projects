#!/bin/bash
#SBATCH --job-name="salixskmer"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --array=1-28 
#SBATCH --mem=32G #this is per job

#--------------------------------------------------------------------------------------------------------------------
#use srun after navigating to scripts directory on /scratch
#run this after activating the env 'skmer'


function subsample_seqtk {

INDIR=$1
OUTDIR=$2

seqtk sample ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz -s 10 20000000 > ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq
seqtk sample ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz -s 10 20000000 > ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq

}

function bbmerge {

INDIR=$1
OUTDIR=$2

#skmer takes uncompressed inputs

merged=${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_merged.fq
unmerged1=${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_unmerged1.fq
unmerged2=${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_unmerged2.fq
#used the hacking method to bbtools into the env directory
$bbtools/bbmerge.sh t=4 in1=${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq in2=${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq out=$merged outu1=$unmerged1 outu2=$unmerged2

#cat all merged and unmerged reads together.  
cat $merged $unmerged1 $unmerged2 > ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.fq 

}




function main {

env_name=skmer
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_salix_28
prefix=salix
####
#constants
easy353path=$HOME/projects/rbge/$USER/env/easy353/bin/Easy353 #this way you don't need to config in .bashrc and thus contain all the packages/path within the env
trimalpath=$HOME/projects/rbge/$USER/env/easy353/bin/trimal/source
bbtools=$HOME/projects/rbge/$USER/env/skmer/bin/bbmap

#subdirectories
ref=$SCRATCH/$WORKDIR/results/353_ref_Salix
RENAMED=$SCRATCH/$WORKDIR/results/renamed
RESULT1=$SCRATCH/$WORKDIR/results/qc1_trimmomatic #this step is omitted eventually
RESULT2=$SCRATCH/$WORKDIR/results/qc2_fastp
RESULT3=$SCRATCH/$WORKDIR/results/phy1_plastome
RESULT4=$SCRATCH/$WORKDIR/results/phy2_alignment_tree
RESULT5=$SCRATCH/$WORKDIR/results/phy3_easy353assembly
RESULT6=$SCRATCH/$WORKDIR/results/phy4_easy353merge
RESULT7=$SCRATCH/$WORKDIR/results/phy5_easy353correctorient
RESULT8=$SCRATCH/$WORKDIR/results/phy6_easy353alignments
RESULT9=$SCRATCH/$WORKDIR/results/phy7_easy353trimal
RESULT10=$SCRATCH/$WORKDIR/results/phy8_easy353concatenate
RESULT11=$SCRATCH/$WORKDIR/results/phy9_easy353phylogeny
RESULT12=$SCRATCH/$WORKDIR/results/phy10_easy353_singlegene_phy
RESULT13=$SCRATCH/$WORKDIR/results/phy11_easy353_high_resolution_phy
RESULT14=$SCRATCH/$WORKDIR/results/phy12_captus_assembly
RESULT15=$SCRATCH/$WORKDIR/results/phy13_captus_extract
RESULT16=$SCRATCH/$WORKDIR/results/phy14_subsampling
RESULT17=$SCRATCH/$WORKDIR/results/phy15_merge_subsample

#USAGES
#=============================================
#Step1: subsampling
subsample_seqtk $RESULT2 $RESULT16
#Step 2 Merge overlapping read pairs**  
bbmerge $RESULT16 $RESULT17



}

main


