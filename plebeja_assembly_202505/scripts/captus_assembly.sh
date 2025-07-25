#!/bin/bash
#SBATCH --job-name="Bcassembly" 
#SBATCH --export=ALL
#SBATCH --partition=medium
#SBATCH --mem=198G 

#
function captus_assembly {

INDIR=$1

cd $INDIR/.. #got to the directory just above the result directory in use
captus clean -r $INDIR #this will generate RESULT2
captus assemble -r $WORKDIR/results/01_clean_reads #the general structure of captus output:RESULT3
cd - #go back to the script directory
}

function captus_extract {

INDIR=$1
REF=$2

cd $INDIR/..
captus extract -a $INDIR -n $REF
}

function main {

#conda create -n captus -c bioconda -c conda-forge captus iqtree -y# T

#UNIVERSAL VARIABLES
env_name=captus 
USER=zedchen 
WORKDIR=${SCRATCH}/plebeja_assembly_202505
prefix=Bp

#CONSTANTS
REF=$HOME/projects/rbge/Begonia_genomes/Reference_Assembelies/conch_genome_v4.fasta

#subdirectories
DATA=$WORKDIR/inputs
RESULT1=$WORKDIR/results/00_qctrimm 
RESULT2=$WORKDIR/results/01_clean_reads #created by captus
RESULT3=$WORKDIR/results/02_assemblies #created by captus
RESULT3=$WORKDIR/results/03_extractions
#$WORKDIR/results/
mkdir $RESULT3

#USAGES
#=============================================
captus_assembly $RESULT1
captus_extract $RESULT3 $REF #reference based assembly

}

main


