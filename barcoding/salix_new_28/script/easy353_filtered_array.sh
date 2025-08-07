#!/bin/bash
#SBATCH --job-name="salix_plastome_phylogeny"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=16
#SBATCH --array=0,10,20,30,40
#SBATCH --mem=32G #change it to the amount of memory you need



function easy353_trimal {

INDIR=$1
OUTDIR=$2

$trimalpath/trimal -in $INDIR/alignment_filter_${SLURM_ARRAY_TASK_ID}.fasta -out $OUTDIR/alignment_trimal_${SLURM_ARRAY_TASK_ID}.fasta -automated1

}

#the phylogeny is resolved poorly???
#the alignment all looks quite accurate?
#wrong sample???
function easy353_phylogeny {

INDIR=$1
OUTDIR=$2

mkdir $OUTDIR/alignment_filter_${SLURM_ARRAY_TASK_ID}
mv $INDIR/alignment_trimal_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/alignment_filter_${SLURM_ARRAY_TASK_ID}
#mafft --maxiterate 10000 $OUTDIR/alignment_filter_${SLURM_ARRAY_TASK_ID}/alignment_trimal_${SLURM_ARRAY_TASK_ID}.fasta > $OUTDIR/alignment_filter_${SLURM_ARRAY_TASK_ID}/alignment_mafft_${SLURM_ARRAY_TASK_ID}.fasta
iqtree -s $OUTDIR/alignment_filter_${SLURM_ARRAY_TASK_ID}/alignment_trimal_${SLURM_ARRAY_TASK_ID}.fasta -bb 10000 -redo -safe -ntmax 16
}


function main {

#UNIVERSAL VARIABLES
env_name=easy353
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_salix_28
prefix=salix
#SET UP DIRECTORIES
easy353path=$HOME/projects/rbge/$USER/env/$env_name/bin/Easy353
trimalpath=$HOME/projects/rbge/$USER/env/$env_name/bin/trimal/source
PARENT=$HOME/scratch/$WORKDIR/results
ref353=$PARENT/353_ref_Salix
#subdirectories
RENAMED=$HOME/scratch/$WORKDIR/results/renamed
CSV=$RENAMED/renamed.csv
RESULT1=$HOME/scratch/$WORKDIR/results/qc1_trimmomatic #this step is omitted eventually
RESULT2=$HOME/scratch/$WORKDIR/results/qc2_fastp
RESULT3=$HOME/scratch/$WORKDIR/results/phy1_plastome
RESULT4=$HOME/scratch/$WORKDIR/results/phy2_alignment_tree
RESULT5=$HOME/scratch/$WORKDIR/results/phy3_easy353assembly
RESULT6=$HOME/scratch/$WORKDIR/results/phy4_easy353merge
RESULT7=$HOME/scratch/$WORKDIR/results/phy5_easy353correctorient
RESULT8=$HOME/scratch/$WORKDIR/results/phy6_easy353alignments
RESULT9=$HOME/scratch/$WORKDIR/results/phy7_easy353trimal
RESULT10=$HOME/scratch/$WORKDIR/results/phy8_easy353concatenate
RESULT11=$HOME/scratch/$WORKDIR/results/phy9_easy353phylogeny
#IN ARRAY

easy353_trimal $RESULT10 $RESULT10  #also rename the contigs to species name
easy353_phylogeny $RESULT10 $RESULT11
}


#execution
main