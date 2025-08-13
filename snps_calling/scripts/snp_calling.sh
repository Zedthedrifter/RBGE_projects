#!/bin/bash
#SBATCH --job-name="Bcassembly" 
#SBATCH --export=ALL
#SBATCH --partition=medium
#SBATCH --mem=3G 
#SBATCH --cpus-per-task=8


function setup_evn {

conda create -n $ENV -y
conda install -c bioconda samtools -y --name $ENV
conda install -c bioconda bcftools -y --name $ENV
conda install -c bioconda vcftools -y --name $ENV
conda install -c bioconda pandas -y --name $ENV #FOR THE EXTRACTION STAGE
}

function setup_workdir {

mkdir $WORKDIR
mkdir $WORKDIR/results
mkdir $INPUTS
mkdir $RESULT1
mkdir $RESULT2

}

#STEP1
function make_vcf {

OUTDIR=$1
INDIR=$2 #THE INPUTDIRECTORY WHERE THE BAM FILES ARE SAVED
ref=$3

# Prep: Converting bam files to VCF files which can then be used as input for the NucBarcoder pipeline
#COPY THE REFERENCE FASTA
rm $INPUTS/REFERENCE.fa* #YOU CAN COMMENT THIS OUT AFTER THE FIRST TIME
cp $ref $INPUTS/REFERENCE.fa #YOU CAN COMMENT THIS OUT AFTER THE FIRST TIME

ls $INDIR/*bam > bamlist.txt #switch this back on later. now i just want to keep three samples:Araucaria-laubenfelsii_NC01610MtM.bam, Araucaria-montana_NC01118P.bam, Araucaria-rulei_NC01190P.bam

bcftools mpileup -Ou -f $INPUTS/REFERENCE.fa --bam-list bamlist.txt | \
    bcftools call -Ou -mv | \
    bcftools filter -s LowQual -e 'QUAL<20 || DP>100' > $OUTDIR/var.flt1.vcf

mv bamlist.txt $RESULT1/bamlist.csv 
#to make the ID_to_scientific_name.csv file
#copy bamlist.csv as ID_to_scientific_name.csv and add the species name manually 
}

#STEP2
function filter_vcf {

INDIR=$1
OUTDIR=$2

# Prep: Filter SNPs using vcftools

vcftools --vcf $INDIR/var.flt1.vcf --out $OUTDIR/var.flt2 --recode --recode-INFO-all \
      --minQ 30 --max-missing 0.95 --minDP 7 --maxDP 100 \
      --min-alleles 2 --max-alleles 2 --remove-indels  --hwe 0.05
      #I removed --mac 5, which requires allele count to be >=5. Allele count is simply the number of times that allele appears over all individuals at that site. this varies for sample size etc. 
}

#STEP3
#YOU NEED TO MAKE A CSV FILE FOR SAMPLE NAMES (JUST COPY THE bamlist.txt FILE IN YOUR CURRENT DIRECTORY AND $RESULT1) AND SPECIES IN $RESULT1, NAME IT ID_to_scientific_name.csv
function extract_loci {

INDIR=$1
OUTDIR=$2
high=$3
low=$4

./calculate_snp.freq.py \
        -v $INDIR/var.flt2.recode.vcf \
        -n $INDIR/ID_to_scientific_name.csv \
        -i $high \
        -l $low \
        -o $OUTDIR #output directory

}
 
#===================================================================

function main {

#parameters and paths to adjust##########################################################################################################################################

#UNIVERSAL VARIABLES
ENV=snps 
WORKDIR=$SCRATCH/SNPS #set up in project, not for anything too big
high=90 #LOWEST FREQUENCY OF THE 'PRESENT' SNP
low=10 #HIGHEST FREQEUENCY OF THE 'ABSENT' SNP


#CONSTANTS
REF=/mnt/shared/projects/rbge/A_projects_Markus/Araucaria/Lib2_mydata_extracted/Araucaria_input-seq_with400Ns_beginend.fas #the fasta file for mapping
data=/mnt/shared/projects/rbge/A_projects_Markus/Araucaria/Lib2_mydata_extracted/exons/21mapped_bwa #directory where the bam files are

#End of parameters and paths to adjust#####################################################################################################################################

#subdirectories
INPUTS=$WORKDIR/results/00_inputs #JUST A PLACE TO PUT COPIES OF REF FILE AND THE FAI INDEXING
RESULT1=$WORKDIR/results/01_raw_data #WHERE THE VCF FILES GO
RESULT2=$WORKDIR/results/02_output #THE SNP STATS AND SELECTED SNPS


#RUN COMMANDS
#SET UP YOUR ENV AND WORK DIRECTORIES
#setup_evn #IF YOU HAVE AN ENV WITH ALL THE REQUIRED TOOLS INSTALLED (SEE THE FUNCTION), YOU CAN SKIP THIS AND ACTIVATE THE CORRESPONDING ENV
#setup_workdir

#CONDA ACTIVATE snps!!!!!!
#conda activate snps BEFORE PROCEEDING!!! THAT'S WHY I COMMENTED OUT THE REST OF THE CODES! 
#STEP1
#make_vcf $RESULT1 $data $REF
#STEP 2
#filter_vcf $RESULT1 $RESULT1
#STEP 3: THIS IS FAST AND YOU DON'T EVEN NEED TO SBATCH IT
extract_loci $RESULT1 $RESULT2 $high $low
}

main


