#!/bin/bash
#SBATCH --job-name="sanger47"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --partition=short
#SBATCH --array=1-352 
#SBATCH --mem=3G #change it to the amount of memory you need

function make_work_dir { 

#mkdir $PARENT

#Create folders under result folder for each program
#QC results


mkdir $RENAMED
mkdir $RESULT1
mkdir $RESULT2
mkdir $RESULT3

}

function rename_contigs {

INDIR=$1
OUTDIR=$2
CSV=$3

./rename_files.py rename_contig -i $INDIR --infile gene_${SLURM_ARRAY_TASK_ID}.fasta --fcsv $CSV -o $OUTDIR
}



function easy353_mafft {

INDIR=$1
OUTDIR=$2

mafft --maxiterate 10000 $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta > $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta
}

function easy353_trimal {

INDIR=$1
OUTDIR=$2
CSV=$3

$trimalpath/trimal -in $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta -out $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta -automated1
./rename_files.py rename_contig -i $OUTDIR --infile gene_${SLURM_ARRAY_TASK_ID}.fasta --fcsv $CSV -o $OUTDIR
}

#make phylogeny for each of the 353 genes
#first step for ASTRAL method
function iqtree_per_gene {

INDIR=$1
OUTDIR=$2

mkdir $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
cp $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
iqtree -s $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta -bb 1000 -redo -safe
}


function main {

#UNIVERSAL VARIABLES
env_name=easy353
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_sanger_47
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
RESULT12=$HOME/scratch/$WORKDIR/results/phy10_easy353_singlegene_phy

#IN ARRAY

#easy353_mafft $RESULT7 $RESULT8
rename_contigs $RESULT8 $RESULT8 $CSV
#easy353_trimal $RESULT8 $RESULT9 $CSV #also rename the contigs to species name
#iqtree_per_gene $RESULT9 $RESULT12
}


#execution
main