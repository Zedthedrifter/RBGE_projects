#!/bin/bash
#SBATCH --job-name="salix47easy353"
#SBATCH --export=ALL
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G 

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

function run_astral {

INDIR=$1
OUTDIR=$2
cvg=$3 #can resolve >= xx phylogeny

echo 'running astral'
./count_mono.py gene_select -i $INDIR -o $OUTDIR --cvg $cvg
#create astral phylogeny with selected genes
java -jar $astral -i $OUTDIR/rsl_${cvg}genes.in.treefile -o $OUTDIR/rsl_${cvg}genes.out.treefile 2>out.log

}

function main {

env_name=easy353
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_sanger_47
prefix=salix
####
easy353path=$HOME/projects/rbge/$USER/env/easy353/bin/Easy353 #this way you don't need to config in .bashrc and thus contain all the packages/path within the env
trimalpath=$HOME/projects/rbge/$USER/env/easy353/bin/trimal/source
astral=$HOME/projects/rbge/$USER/env/easy353/bin/ASTRAL/astral.5.7.8.jar #the java script that build an astral tree from a collection of trees
PARENT=$HOME/scratch/$WORKDIR/results
ref=$PARENT/353_ref_Salix
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
#USAGES
#=============================================
#make work directories

#mkdir $RESULT5
#mkdir $RESULT6
#mkdir $RESULT7
#mkdir $RESULT8
#mkdir $RESULT9
#mkdir $RESULT10
#mkdir $RESULT11
#mkdir $RESULT13


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
#./count_mono.py process_treefiles -i $RESULT12 -o $RESULT13

#use astral method to make phylogeny: how many high-resolution genes per taxa do we need?
run_astral $RESULT12 $RESULT13 2 #get at least 2 genes for each mono taxa
run_astral $RESULT12 $RESULT13 3 #get at least 3 genes for each mono taxa
run_astral $RESULT12 $RESULT13 4 #get at least 4 genes for each mono taxa
run_astral $RESULT12 $RESULT13 5 #get at least 5 genes for each mono taxa
run_astral $RESULT12 $RESULT13 6 #get at least 6 genes for each mono taxa

echo $genes
}

main

