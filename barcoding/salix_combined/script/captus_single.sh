#!/bin/bash
#SBATCH --job-name="slxcaptus"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --partition=long
#SBATCH --mem=50G #change it to the amount of memory you need

function captus_assembly {

INDIR=$1

cd $INDIR/.. #got to the directory just above the result directory in use
captus clean -r $INDIR #this will generate RESULT2 #this part is done
captus assemble -r $WORKDIR/results/01_clean_reads #the general structure of captus output:RESULT3
cd - #go back to the script directory
}

function main {

#UNIVERSAL VARIABLES
env_name=easy353
username=zedchen 
WORKDIR=Barcoding_km/barcoding_combined_salix #the last time i just ran it with salix_47 by mistake
prefix=salix
reads=$HOME/projects/rbge/$username/barcoding/salix_combined/results/00_reads #better keeping the raw reads on project
####
easy353path=$HOME/projects/rbge/$username/env/easy353/bin/Easy353 #this way you don't need to config in .bashrc and thus contain all the packages/path within the env
trimalpath=$HOME/projects/rbge/$username/env/easy353/bin/trimal/source
astral=$HOME/projects/rbge/$username/env/easy353/bin/ASTRAL/astral.5.7.8.jar #the java script that build an astral tree from a collection of trees
PARENT=$HOME/scratch/$WORKDIR/results
ref=$PARENT/353_ref_Salix
#subdirectories
RESULT0=$reads
RESULT1=$SCRATCH/$WORKDIR/results/01_easy353_compile
RESULT2=$SCRATCH/$WORKDIR/results/02_easy353_alignment
RESULT3=$SCRATCH/$WORKDIR/results/03_easy353_singlegene_phy
RESULT4=$SCRATCH/$WORKDIR/results/04_easy353_astral
#captus
RESULT5=$SCRATCH/$WORKDIR/results/01_clean_reads
#RUN FUNCTIONS 
captus_assembly $RESULT0

}


#execution
main