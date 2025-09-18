#!/bin/bash
#SBATCH --job-name="easy353"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --array=1
#SBATCH --mem=50G #change it to the amount of memory you need

function rename_cp { #a much quicker way to do it

INDIR=$1
OUTDIR=$2
prefix=$3

cp $INDIR/45926\#${SLURM_ARRAY_TASK_ID}_R1.f*gz $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz
cp $INDIR/45926\#${SLURM_ARRAY_TASK_ID}_R2.f*gz $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz

}

function fastp_qc {

#SKIP TRIM ADAPTOR
INDIR=$1
OUTDIR=$2

fastp -i ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz -I ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
       -o ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_1.fq.gz -O ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_2.fq.gz \
       --unpaired1 ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_unpaired1.fq.gz --unpaired2 ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_unpaired2.fq.gz \
       -j ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.json \
       -h ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.html \
       --dont_overwrite --detect_adapter_for_pe -c > ${OUTDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}.fastp.log


}

function easy353_assembly {

INDIR=$1
OUTDIR=$2
ref353=$3

$easy353path/easy353.py -1 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz -2 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz \
                        -r $ref353/353gene -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID} \
                        -fk 31 -ak 41 -ft 1 -at 4 -kmer_limit 8 


}

#build astral trees 
function run_astral {

INDIR=$1
OUTDIR=$2
cvg=$3 #can resolve >= xx phylogeny

echo 'running astral'
./count_mono.py gene_select -i $INDIR -o $OUTDIR --cvg $cvg
#create astral phylogeny with selected genes
java -jar $astral -i $OUTDIR/rsl_${cvg}genes.in.treefile -o $OUTDIR/rsl_${cvg}genes.out.treefile 2>out.log

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


#make phylogeny for each of the 353 genes
#first step for ASTRAL method
function iqtree_per_gene {

INDIR=$1
OUTDIR=$2

mkdir $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
cp $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
#remove bad taxa (sample is bad) 'repens' and other all gap contigs
./remove_taxa.py $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta Not_removing_taxa 
#remove the outliers from an alignment
#remove sequences less than 65% identical
./remove_highlyhetero.py $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta -t 0.65 -m majority 
iqtree -s $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta -bb 1000 -redo -safe --mem 8G
}


function make_work_dir { 

mkdir $WORKDIR
mkdir $RESULT1
mkdir $RESULT2
mkdir $RESULT3
mkdir $RESULT4
mkdir $RESULT5
mkdir $RESULT6
mkdir $RESULT7
}

function main {

#UNIVERSAL VARIABLES
env_name=easy353
USER=zedchen 
WORKDIR=$SCRATCH/Barcoding_km/barcoding_combined_salix/results
DATA=/mnt/shared/projects/rbge/pholling/barcoding/salix_combined/results/00_reads
CSV=/mnt/shared/projects/rbge/pholling/barcoding/salix_combined/results/00_reads/renamed.csv
prefix=salix
#SET UP DIRECTORIES
easy353path=$HOME/apps/manual/Easy353
astral=$HOME/apps/manual/ASTRAL/astral.5.7.8.jar #the java script that build an astral tree from a collection of trees
ref353=/mnt/shared/scratch/zchen/Barcoding_km/barcoding_sanger_47/results/353_ref_Salix
#subdirectories

RESULT1=$WORKDIR/01_easy353assembly #this step is omitted eventually
RESULT2=$WORKDIR/02_easy_merge
RESULT3=$WORKDIR/03_easy_reorient
RESULT4=$WORKDIR/04_easy_rename
RESULT5=$WORKDIR/05_easy_aligned
RESULT6=$WORKDIR/06_easy_treefiles
RESULT7=$WORKDIR/07_easy_phylogeny
RESULT=$WORKDIR/0
#IN ARRAY
make_work_dir

function single1 {
#merge the recovered genes from different samples
cd $RESULT1 #the design of the code means that you must be in the result directory for the files to be found
python $easy353path/script/multi_sample_gene_merge.py -i $RESULT1/* -p 0 -o $RESULT2
#Correct Orientation
python $easy353path/script/correct_seq_ori.py -i $RESULT2/combine -r $ref353/353gene -o $RESULT3
cd -
#rename the contigs in fasta
./rename_files.py rename_fasta_easy353 -i $RESULT3 -o $RESULT3
}

function single2 {

#use astral method to make phylogeny: how many high-resolution genes per taxa do we need?
run_astral $RESULT6 $RESULT7 2 #get at least 2 genes for each mono taxa
run_astral $RESULT6 $RESULT7 3 #get at least 3 genes for each mono taxa
run_astral $RESULT6 $RESULT7 4 #get at least 4 genes for each mono taxa
run_astral $RESULT6 $RESULT7 5 #get at least 5 genes for each mono taxa
run_astral $RESULT6 $RESULT7 6 #get at least 6 genes for each mono taxa
}

#rename_cp $DATA $RENAMED $prefix
#fastp_qc $RENAMED $RESULT2 #

#STEP_1: ASSEMBLE: 1-86%12
#you still need to assemble all of them: the renamed genes might not compile very well together
#easy353_assembly $DATA $RESULT1 $ref353 #use fastp cleaned reads

#STEP 2: MERGE ASSEMBLY, REORIENT, RENAME GENE FILES
#RUN AS SINGLE: 1

#STEP 3 SINGLE: REORIENT AND CORRECT,
#single1

#STEP 4: IQTREE: GENE ARRAY: 1-352
#rename_contigs $RESULT3 $RESULT4 $CSV
#easy353_mafft $RESULT4 $RESULT5
#iqtree_per_gene $RESULT5 $RESULT6
#

#STEP5 ASTRAL
#after running phylogeny of individual genes, check how many monophyletic species are resolved and select genes of high resolution
#./count_mono.py process_treefiles -i $RESULT6 -o $RESULT7

single2
}


#execution
main