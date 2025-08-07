#!/bin/bash
#SBATCH --job-name="slxastral"
#SBATCH --export=ALL
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G 

#UNIVERSAL VARIABLES
env_name=easy353
USER=zedchen 
WORKDIR=Barcoding_km # subdirectory on scratch where the project will be carried out: $HOME/scratch/$work_dir/

#conda create -n $env_name #make an env 
#--------------------------------------------------------------------------------------------------------------------
#use srun after navigating to scripts directory on /scratch
#run this after activating the env 'salix'


#Install samtools, NanoPlot, porechop, chopper
#these are the old QC packages from amplicon-seq project, some might be useful, such as samtools, biopython,etc
function install_basic {    

#install to specific env
conda install -c bioconda samtools -y --name $env_name
conda install pip -y --name $env_name #install pip
}

function pip_install {

pip=$HOME/projects/rbge/$USER/env/$env_name/bin/pip #install with pip
#$pip install psutil 
#$pip install biopython
#$pip install matplotlib
#$pip install requests 
#$pip install beautifulsoup4
#$pip install amas
$pip install seaborn 
}

#--------------------------------------------------------------------------------------------------------------------
#these are packages for the salix project

function install_easy353 {

git clone https://github.com/plant720/Easy353.git $HOME/projects/rbge/$USER/env/$env_name/bin/Easy353
cd $HOME/projects/rbge/$USER/env/$env_name/bin/Easy353
python3 setup.py install --user
}

function install_trimal { 

mkdir $HOME/projects/rbge/$USER/env/$env_name/bin/trimal
git clone https://github.com/inab/trimal.git $HOME/projects/rbge/$USER/env/$env_name/bin/trimal
cd $HOME/projects/rbge/$USER/env/$env_name/bin/trimal/source
make

}

function install_tree {

conda install -c bioconda iqtree -y --name $env_name
conda install -c bioconda figtree -y --name $env_name
conda install -c conda-forge mafft -y --name $env_name
conda install -c bioconda orthofinder -y --name $env_name
conda install -c anaconda scipy -y --name $env_name
}

function build_353_db {

outdir=$1
genus=$2

#config and add permission
chmod 700 $easy353path/*py
#build database
$easy353path/build_database.py -o $outdir -c $genus -t 10 -generate  

}


#make directories for storing results
#scripts and data directories are already made

function setup_env {

mkdir $easy353path
install_easy353
install_basic
install_tree
pip_install
install_trimal #trim alignment: remove regions of poor alignment: https://vicfero.github.io/trimal/index.html#installation_sec
}

function run_astral_coverage {

INDIR=$1
OUTDIR=$2
cvg=$3 #can resolve >= xx phylogeny

echo 'running astral'
./count_mono.py gene_select -i $INDIR -o $OUTDIR --cvg $cvg
#create astral phylogeny with selected genes
java -jar $astral -i $OUTDIR/rsl_${cvg}genes.in.treefile -o $OUTDIR/rsl_${cvg}genes.out.treefile 2>out.log

}

function run_astral_manual {

INDIR=$1
OUTDIR=$2
genelist=$3 #sets of genes, combinations

echo 'running astral with manually selected genes'
./count_mono.py gene_manual -i $INDIR -o $OUTDIR --genes $genelist
#create astral phylogeny with selected genes
for infile in $(ls $OUTDIR/manual_*.in.treefile) #check all the manual trees
do
java -jar $astral -i $infile -o ${infile/in.treefile/out.treefile} 2>out.log
done
}

function main {

env_name=easy353
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_combined_salix #the last time i just ran it with salix_47 by mistake
prefix=salix
####
easy353path=$HOME/projects/rbge/$USER/env/easy353/bin/Easy353 #this way you don't need to config in .bashrc and thus contain all the packages/path within the env
trimalpath=$HOME/projects/rbge/$USER/env/easy353/bin/trimal/source
astral=$HOME/projects/rbge/$USER/env/easy353/bin/ASTRAL/astral.5.7.8.jar #the java script that build an astral tree from a collection of trees
PARENT=$HOME/scratch/$WORKDIR/results
ref=$PARENT/353_ref_Salix
#subdirectories
RESULT1=$SCRATCH/$WORKDIR/results/01_easy353_compile
RESULT2=$SCRATCH/$WORKDIR/results/02_easy353_alignment
RESULT3=$SCRATCH/$WORKDIR/results/03_easy353_singlegene_phy
RESULT4=$SCRATCH/$WORKDIR/results/04_easy353_astral
#USAGES
#=============================================
#make work directories

#setup_env #the function that actualy set up env, already did with salix, reuse env, no need to run

#build_353_db $ref $prefix
#merge the recovered genes from different samples
#cd $RESULT5 #the design of the code means that you must be in the result directory for the files to be found
#python $easy353path/script/multi_sample_gene_merge.py -i $RESULT5/* -p 0 -o $RESULT6
#Correct Orientation
#python $easy353path/script/correct_seq_ori.py -i $RESULT6/combine -r $ref/353gene -o $RESULT7
#cd $HOME/scratch/$WORKDIR/scripts
#rename the contigs in fasta
#./rename_files.py rename_fasta_easy353 -i $RESULT7 -o $RESULT7
#./rename_files.py gene_alignment_summary -i $RESULT9 -o $RESULT9
#./rename_files.py easy353_salix_filter -i $RESULT9 -o $RESULT10

#after running phylogeny of individual genes, check how many monophyletic species are resolved and select genes of high resolution
#./count_mono.py process_treefiles -i $RESULT3 -o $RESULT4

#use astral method to make phylogeny: how many high-resolution genes per taxa do we need?
#run_astral_coverage $RESULT3 $RESULT4 2 #get at least 2 genes for each mono taxa
#run_astral_coverage $RESULT3 $RESULT4 3 #get at least 3 genes for each mono taxa
#run_astral_coverage $RESULT3 $RESULT4 4 #get at least 4 genes for each mono taxa
#run_astral_coverage $RESULT3 $RESULT4 5 #get at least 5 genes for each mono taxa
#run_astral_coverage $RESULT3 $RESULT4 6 #get at least 6 genes for each mono taxa

#manually select a set of genes
run_astral_manual $RESULT3 $RESULT4 $RESULT4/manual_gene_sets.txt 

echo $genes
}

main
