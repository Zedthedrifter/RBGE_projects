#!/bin/bash
#SBATCH --job-name="salixskmer"
#SBATCH --export=ALL
#SBATCH --partition=medium
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
conda create -n $env_name #make an env 
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

function install_bbtools {

#wget https://downloads.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz --directory-prefix=$HOME/projects/rbge/$USER/env/$env_name/bin
#cd $HOME/projects/rbge/$USER/env/$env_name/bin
#tar -xvzf $HOME/projects/rbge/$USER/env/$env_name/bin/BBMap_39.01.tar.gz
rm $HOME/projects/rbge/$USER/env/$env_name/bin/BBMap_39.01.tar.gz
}


function install_tree {

conda install -c bioconda iqtree -y --name $env_name
conda install -c bioconda figtree -y --name $env_name
conda install -c conda-forge mafft -y --name $env_name
conda install -c bioconda orthofinder -y --name $env_name
conda install -c anaconda scipy -y --name $env_name
}



#make directories for storing results
#scripts and data directories are already made

function setup_env {

#conda create -n $env_name
#conda install -c bioconda samtools -y --name $env_name
#conda install -c bioconda seqtk -y --name $env_name
#conda install pip -y --name $env_name #install pip
#conda install -c bioconda skmer=3.3.0 -y --name $env_name
conda install -c bioconda fastme=2.1.6.3 -y --name $env_name
#install_bbtools
#some installing functions
#install_tree
#pip_install
}

function skmer_ref {

INDIR=$1
OUTDIR=$2

skmer reference $INDIR -s 100000 -S 42 -p 24 -t -o $OUTDIR

}

function skmer_subsample {

INDIR=$1
OUTDIR=$2

#subsampling to make 100 replicates
skmer subsample -b 100 $INDIR -s 100000 -S 42 -p 24 -t -i 0 -sub $OUTDIR
}

function skmer_correct {

INDIR1=$1
INDIR2=$2
#Correct estimates:  
skmer correct -main $INDIR1/dimtrx_main.txt -sub $INDIR2
}


function main {

env_name=skmer
USER=zedchen 
WORKDIR=Barcoding_km/barcoding_salix_28
prefix=salix
####
#constants
easy353path=$HOME/projects/rbge/$USER/env/easy353/bin/Easy353 #this way you don't need to config in .bashrc and thus contain all the packages/path within the env
trimalpath=$HOME/projects/rbge/$USER/env/easy353/bin/trimal/source
bbtools=$HOME/projects/rbge/$USER/env/skmer/bin/bbmap

#subdirectories
ref=$SCRATCH/$WORKDIR/results/353_ref_Salix
RENAMED=$SCRATCH/$WORKDIR/results/renamed
RESULT1=$SCRATCH/$WORKDIR/results/qc1_trimmomatic #this step is omitted eventually
RESULT2=$SCRATCH/$WORKDIR/results/qc2_fastp
RESULT3=$SCRATCH/$WORKDIR/results/phy1_plastome
RESULT4=$SCRATCH/$WORKDIR/results/phy2_alignment_tree
RESULT5=$SCRATCH/$WORKDIR/results/phy3_easy353assembly
RESULT6=$SCRATCH/$WORKDIR/results/phy4_easy353merge
RESULT7=$SCRATCH/$WORKDIR/results/phy5_easy353correctorient
RESULT8=$SCRATCH/$WORKDIR/results/phy6_easy353alignments
RESULT9=$SCRATCH/$WORKDIR/results/phy7_easy353trimal
RESULT10=$SCRATCH/$WORKDIR/results/phy8_easy353concatenate
RESULT11=$SCRATCH/$WORKDIR/results/phy9_easy353phylogeny
RESULT12=$SCRATCH/$WORKDIR/results/phy10_easy353_singlegene_phy
RESULT13=$SCRATCH/$WORKDIR/results/phy11_easy353_high_resolution_phy
RESULT14=$SCRATCH/$WORKDIR/results/phy12_captus_assembly
RESULT15=$SCRATCH/$WORKDIR/results/phy13_captus_extract
RESULT16=$SCRATCH/$WORKDIR/results/phy14_subsample
RESULT17=$SCRATCH/$WORKDIR/results/phy15_merge_subsample
RESULT18=$SCRATCH/$WORKDIR/results/phy16_dimtrx_main #skmer distance matrix
RESULT19=$SCRATCH/$WORKDIR/results/phy17_skmersubsample
#USAGES
#=============================================
#setup_env #the function that actualy set up env
skmer_ref $RESULT17 $RESULT18 
skmer_subsample $RESULT17 $RESULT19
skmer_correct $RESULT18 $RESULT19
}

main


