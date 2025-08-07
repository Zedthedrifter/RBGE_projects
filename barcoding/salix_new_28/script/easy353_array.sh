#!/bin/bash
#SBATCH --job-name="salix_plastome_phylogeny"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --array=1-28 
#SBATCH --mem=32G #this is per job

function make_work_dir { 

#mkdir $PARENT

#Create folders under result folder for each program
#QC results


mkdir $RENAMED
mkdir $RESULT1
mkdir $RESULT2
mkdir $RESULT3

}


#Whole plastome assembly: array
function plastome_assembly {

INDIR=$1
OUTDIR=$2
ref=$3

get_organelle_from_reads.py -s $ref -1 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz -2 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
                            -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID} -R 15 -k 21,45,65,85,105 -F embplant_pt

}

function easy353_assembly {

INDIR=$1
OUTDIR=$2
ref353=$3

$easy353path/easy353.py -1 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz -2 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
                        -r $ref353/353gene -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID} \
                        -fk 31 -ak 41 -ft 1 -at 4 -kmer_limit 8 


}

function captus_assembly {

INDIR=$1
OUTDIR1=$2
OUTDIR2=$3
#i might need to break this into multiple smaller functions? already installed into easy353 env
cd $PARENT
captus assemble -r ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz -o $OUTDIR1 #could have submitted as a directly, but this way it can run easily in parallel
                
captus extract -a ${OUTDIR1}/${prefix}_${SLURM_ARRAY_TASK_ID}__captus-asm -n Angiosperms353 -o $OUTDIR2               
}


function main {

#UNIVERSAL VARIABLES
env_name=easy353
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_salix_28
prefix=salix
#SET UP DIRECTORIES
easy353path=$HOME/projects/rbge/$USER/env/$env_name/bin/Easy353
PARENT=$HOME/scratch/$WORKDIR/results
ref353=$PARENT/353_ref_Salix
#subdirectories
RENAMED=$HOME/scratch/$WORKDIR/results/renamed
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
RESULT13=$HOME/scratch/$WORKDIR/results/phy11_easy353_high_resolution_phy
RESULT14=$HOME/scratch/$WORKDIR/results/phy12_captus_assembly
RESULT15=$HOME/scratch/$WORKDIR/results/phy13_captus_extract
#IN ARRAY
#plastome_assembly $RESULT2 $RESULT3 $ref #you have to write the inputs
#easy353_assembly $RESULT2 $RESULT5 $ref353 #use fastp cleaned reads
captus_assembly $RESULT2 $RESULT14 $RESULT15
}


#execution
main