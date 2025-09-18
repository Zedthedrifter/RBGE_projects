#!/bin/bash
#SBATCH --job-name="potabusco"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --mem=10G 

#but a few files actually doesn't exist



#make directories
function setup_dir {

mkdir $WORKDIR
mkdir $WORKDIR/results
mkdir $WORKDIR/refs
mkdir $REFR1

}

#run MuMmer
function chromosome_homolog {

reference=$1
query=$2
OUTDIR=$3

mkdir $OUTDIR
cd $OUTDIR
#
dnadiff $reference $query

}
function main {

#UNIVERSAL VARIABLES=======================================================================================================================
env_name=snps
USER=zedchen 
WORKDIR=$SCRATCH/Barcoding_km/SNP_potamogeton #actually this is salix. need to rename
prefix=salix
#SLURM_ARRAY_TASK_ID=21

####constants############
#DATA=$HOME/projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/
#==========================================================================================================================================


#subdirectories
REF=$WORKDIR/refs
REFR1=$WORKDIR/results/ref_01_blast #REF RESULT


#USAGES
#=============================================
setup_dir #WILL NOT OVERWRITE DIR, JUST LEAVE IT

#REFERENCE CURATION
#./process_ref.py $REF $REF
#chromosome_homolog: 10G,20G,60G mem
chromosome_homolog $WORKDIR/refs/Sh_chromosome_1.fasta  \
                   $WORKDIR/refs/chromosome_1.fasta \
                   $REFR1/chromosome1 #
 
}

main


