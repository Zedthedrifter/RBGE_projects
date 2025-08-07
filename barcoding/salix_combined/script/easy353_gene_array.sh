#!/bin/bash
#SBATCH --job-name="salixall"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --partition=short
#SBATCH --array=1-352
#SBATCH --mem=2G #change it to the amount of memory you need

function collect_contigs {

OUTDIR=$1

#put the genes from all three folders into the same file, so that now we have 46+28+13 inputs (well... some samples were sequenced twice)
cat $SCRATCH/Barcoding_km/barcoding_sa*/results/phy6_easy353alignments/gene_${SLURM_ARRAY_TASK_ID}.fasta > $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta #collect all the easy353 assembly results

}


function easy353_mafft {

INDIR=$1
OUTDIR=$2

mafft --maxiterate 10000 $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta > $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta
}


#make phylogeny for each of the 353 genes
#first step for ASTRAL method
function iqtree_per_gene {

INDIR=$1
OUTDIR=$2

mkdir $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
cp $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
#remove bad taxa (sample is bad) 'repens' and other all gap contigs
./remove_taxa.py $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta repens 
#remove the outliers from an alignment
#remove sequences less than 65% identical
./remove_highlyhetero.py $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta -t 0.65 -m majority 
iqtree -s $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta -bb 1000 -redo -safe --mem 8G
}



function main {

#UNIVERSAL VARIABLES
env_name=easy353
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_combined_salix
prefix=salix
#SET UP DIRECTORIES
easy353path=$HOME/projects/rbge/$USER/env/$env_name/bin/Easy353
trimalpath=$HOME/projects/rbge/$USER/env/$env_name/bin/trimal/source
PARENT=$SCRATCH/$WORKDIR/results
ref353=$PARENT/353_ref_Salix
#subdirectories
RENAMED=$SCRATCH/$WORKDIR/results/renamed
#CSV=$RENAMED/renamed.csv #this one needs to be remake for later steps
RESULT1=$SCRATCH/$WORKDIR/results/01_easy353_compile
RESULT2=$SCRATCH/$WORKDIR/results/02_easy353_alignment
RESULT3=$SCRATCH/$WORKDIR/results/03_easy353_singlegene_phy
RESULT4=$SCRATCH/$WORKDIR/results/04_easy353_astral

#set up directories

#IN ARRAY
#compile the easy353 assembly results
#collect_contigs $RESULT1
#align the contigs again
#easy353_mafft $RESULT1 $RESULT2
#generate treefiles for individual genes
iqtree_per_gene $RESULT2 $RESULT3
#go to single.sh to process all trees: count_mono.py etc. 
}


#execution
main