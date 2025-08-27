#!/bin/bash
#SBATCH --job-name="Gpneu"
#SBATCH --export=ALL
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G 

function download_data {

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR9124426/ERR9124426 -O $DATA/alpine
$sratools/fasterq-dump $DATA/alpine --split-files --progress --threads 4 -O $DATA
#
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR5554463/ERR5554463 -O $DATA/norway
$sratools/fasterq-dump $DATA/norway --split-files --progress --threads 4 -O $DATA

rm $DATA/alpine $DATA/norway -f
#rename the files for captus assembly
cd $DATA
for i in $(ls *1.fastq); do mv $i ${i/1.fastq/R1.fastq}; done
for i in $(ls *2.fastq); do mv $i ${i/2.fastq/R2.fastq}; done

}

#===========================================================================

#CAPTUS: clean reads and assembly
function captus_assembly {

INDIR=$1

cd $INDIR/.. #got to the directory just above the result directory in use
#captus clean -r $INDIR #
#mem=30G
captus assemble -r $WORKDIR/results/01_clean_reads \
                --k_list 31,39,47,63,79,95 \
                --overwrite \
                --ram 30  

#read length=101. reduce kmer size to 95 at most
#the general structure of captus output:RESULT3
cd - #go back to the script directory
}

#===========================================================================

#BUSCO
function busco_assess {

INFILE=$1
OUTDIR=$2

busco -i $INFILE -o $OUTDIR -m genome -c 32 -l eudicotyledons_odb12 -r
}

#===========================================================================

#QUAST de novo assembly evaluation
function quast_assess {

INFILE=$1
OUTDIR=$2

quast.py $INFILE -o $OUTDIR  
}


#===========================================================================

#SSR detection
function MISA_SSR {

INFILE=$1

#remove reads <400 bp
seqtk seq -L 400 $INFILE > ${INFILE/fasta/filtered.fa}
#
misa.pl ${INFILE/fasta/filtered.fa}
}
#===========================================================================

function setup_dir {

mkdir $DATA
mkdir $WORKDIR/results
mkdir $RESULT3
mkdir $RESULT4
mkdir $RESULT5

}

#===========================================================================

function main {

#conda create -n captus -c bioconda -c conda-forge captus iqtree -y# T

#UNIVERSAL VARIABLES
env_name=captus 
USER=zedchen 
WORKDIR=$SCRATCH/G_pneumonanthe
prefix=Gp
flank_size=150

#CONSTANTS
sratools=$HOME/apps/manual/sratoolkit.3.2.1-ubuntu64/bin/
pip=$HOME/apps/env/easy353/bin/pip
STAR=$HOME/apps/RNA_polish/bin/STAR/bin/Linux_x86_64/STAR
#subdirectories
DATA=$WORKDIR/results/data
#RNA=
RESULT1=$WORKDIR/results/01_clean_reads #created by captus
RESULT2=$WORKDIR/results/02_assemblies #created by captus
RESULT3=$WORKDIR/results/03_BUSCO
RESULT4=$WORKDIR/results/04_SSRs
RESULT5=$WORKDIR/results/05_

setup_dir

#STEP_1 DOWNLOAD
#download_data

#STEP_2 CAPTUS ASSEMBLY
#captus_assembly $DATA
#MANUALLY MOVE
#for i in $(ls $RESULT2|grep captus-asm); do mv $RESULT2/$i/*/assembly.fasta $RESULT4/${i/__captus-asm/}.fasta; done

#STEP_3 BUSCO
#run in env Busco
#busco_assess $RESULT4/..... $RESULT3

#quast_assess $infile $

#STEP_4 SSR detection: run in base. seqtk is installed in base and MISA is universal
#for i in $(ls $RESULT4/*.fasta); do MISA_SSR $i; done #files will be saved to the same directory as the input file

#STEP_5: SEQ EXTRACTION
#flanking size = 100
#for INFILE in $(ls $RESULT2|grep captus-asm); \
#   do ./extract_flanked.py $RESULT4/${INFILE/__captus-asm/.filtered.fa} $RESULT4/${INFILE/__captus-asm/.filtered.fa.misa} ${INFILE/__captus-asm/} $flank_size; \
#   done

function detect_conserved {

flank=$1

STEP_5: IDENTIFY CONSERVED SSR #this one runs ./extract_flanked.py 
./detect_conserved_SSR.py  $RESULT4/alpine.filtered.fa $RESULT4/alpine.filtered.fa.misa alpine \
                           $RESULT4/norway.filtered.fa $RESULT4/norway.filtered.fa.misa norway \
                           $flank
}

detect_conserved 150
#detect_conserved 100

}

main