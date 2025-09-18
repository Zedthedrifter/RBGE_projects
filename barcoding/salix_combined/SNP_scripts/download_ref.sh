#!/bin/bash
#SBATCH --job-name="potabusco"
#SBATCH --export=ALL
#SBATCH --partition=medium
#SBATCH --mem=2G 

#but a few files actually doesn't exist

function download_caprea {

OUTFILE=$1

wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037743.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037744.1 -O - >> $OUTFILE 
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037745.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037746.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037747.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037748.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037749.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037750.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037751.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037752.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037753.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037754.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037755.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037756.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037757.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037758.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037759.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037760.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ037761.1 -O - >> $OUTFILE

}

function download_cinerea {

OUTFILE=$1

wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253417.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253418.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253419.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253420.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253421.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253422.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253423.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253424.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253425.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253426.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253427.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253428.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253429.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253430.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253431.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253432.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253433.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ253434.1 -O - >> $OUTFILE

}

function download_reticulata {

OUTFILE=$1

wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287372.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287373.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287374.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287375.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287376.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287377.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287378.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287379.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287380.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287381.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287382.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287383.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287384.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287385.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287386.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287387.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287388.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287389.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287390.1 -O - >> $OUTFILE

}

function download_viminalis {

OUTFILE=$1

wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287308.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287309.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287310.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287311.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287312.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287313.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287314.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287315.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287316.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287317.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287318.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287319.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287320.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287321.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287322.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287323.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287324.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287325.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ287326.1 -O - >> $OUTFILE

}

function download_herbacea {

OUTFILE=$1

wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286651.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286652.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286653.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286654.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286655.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286656.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286657.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286658.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286659.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286660.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286661.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286662.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286663.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286664.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286665.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286666.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286667.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286668.1 -O - >> $OUTFILE
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OZ286669.1 -O - >> $OUTFILE

}

function download_all_ref {

OUTDIR=$1

#download_caprea $OUTDIR/S_caprea.chromosome.fasta
download_cinerea $OUTDIR/S_cinerea.chromosome.fasta
download_reticulata $OUTDIR/S_reticulata.chromosome.fasta
download_viminalis $OUTDIR/S_viminalis.chromosome.fasta
download_herbacea $OUTDIR/S_herbacea.chromosome.fasta
}

#make directories
function setup_dir {

mkdir $WORKDIR
mkdir $WORKDIR/results
mkdir $WORKDIR/refs


}

function main {

#UNIVERSAL VARIABLES=======================================================================================================================
env_name=busco
USER=zedchen 
WORKDIR=$SCRATCH/Barcoding_km/SNP_potamogeton #
prefix=potamogeton
#SLURM_ARRAY_TASK_ID=21

####constants############
DATA=$HOME/projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/
#==========================================================================================================================================


#subdirectories
REF=$WORKDIR/refs
RESULT1=$WORKDIR/results/01_clean_reads
RESULT2=$WORKDIR/results/02_captus_assembly
RESULT3=$WORKDIR/results/03_captus_extract


#USAGES
#=============================================
setup_dir
download_all_ref $REF


}

main


