#!/bin/bash
#SBATCH --job-name="salixeasy353"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G 



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


function install_tree {

conda install -c bioconda iqtree -y --name $env_name
conda install -c bioconda figtree -y --name $env_name
conda install -c conda-forge mafft -y --name $env_name
conda install -c bioconda orthofinder -y --name $env_name
conda install -c anaconda scipy -y --name $env_name
}



function setup_env {

install_basic
install_tree
#install_easy353
}

function main {

#UNIVERSAL VARIABLES
env_name=busco
USER=zedchen 
WORKDIR=Barcoding_km # subdirectory on scratch where the project will be carried out: $HOME/scratch/$work_dir/
prefix=salix
####
easy353path=$HOME/projects/rbge/$USER/env/$env_name/bin/Easy353 #this way you don't need to config in .bashrc and thus contain all the packages/path within the env
trimalpath=$HOME/projects/rbge/$USER/env/$env_name/bin/trimal/source
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
#USAGES
#=============================================
setup_env #the function that actualy set up env




}

main


