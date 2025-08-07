#!/bin/bash
#SBATCH --job-name="salixiqtree"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G 

function make_work_dir { 

#mkdir $PARENT

#Create folders under result folder for each program
#QC results


mkdir $RESULT4

}


function phylogeny {

INDIR=$1
OUTDIR=$2
CSV=$3

#extract the graph1 contigs from fasta file, named to the same prefix as in slurm array for consistency. need to run only once
#./rename_files.py rename_fasta -i $INDIR -o $OUTDIR --outfile ${prefix}_plastomes.fasta
#./rename_files.py rename_contig -i $OUTDIR --infile ${prefix}_plastomes.fasta --fcsv $CSV -o $OUTDIR
mafft --maxiterate 10000 $OUTDIR/${prefix}_plastomes.fasta > $OUTDIR/${prefix}_plastomes.aligned.fasta
#./rename_files.py rename_contig -i $OUTDIR --infile ${prefix}_plastomes.aligned.fasta --fcsv $CSV -o $OUTDIR
iqtree -s $OUTDIR/${prefix}_plastomes.aligned.fasta -bb 1000 -redo -safe
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

#IN ARRAY
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >> plastome_debug.log
#run mafft and iqtree
phylogeny $RESULT3 $RESULT4 $RENAMED/renamed.csv
}


#execution
main