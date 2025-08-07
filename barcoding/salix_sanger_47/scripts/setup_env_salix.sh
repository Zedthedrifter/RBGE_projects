#!/bin/bash
#KEY UNIVERSAL VARIABLES

env_name=salix
USER=zedchen 
WORKDIR=Barcoding_km # subdirectory on scratch where the project will be carried out: $HOME/scratch/$work_dir/

#conda create -n $env_name #make an env 
#--------------------------------------------------------------------------------------------------------------------
#use srun after navigating to scripts directory on /scratch
#run this after activating the env 'salix'


#Install samtools, NanoPlot, porechop, chopper
#these are the old QC packages from amplicon-seq project, some might be useful, such as samtools, biopython,etc
function install_amp {    

#install to specific env
conda install -c bioconda samtools -y --name $env_name
conda install -c bioconda nanoplot -y --name $env_name
conda install -c bioconda porechop -y --name $env_name
mamba install chopper
mamba update chopper 
#Install dependencies of ampliconsorter 
#python3 -m pip install edlib
#python3 -m pip install biopython
#python3 -m pip install matplotlib
#/bin/python3: No module named pip
#Install dependencies of ampliconsorter 
#
conda install pip -y --name $env_name #install pip
}

function pip_install {

pip=$HOME/projects/rbge/$USER/env/$env_name/bin/pip #install with pip
$pip install edlib
$pip install biopython
$pip install matplotlib
}

#--------------------------------------------------------------------------------------------------------------------
#these are packages for the salix project

function install_QC {

#conda install -c bioconda fastp -y --name $env_name
conda install -c bioconda trimmomatic -y --name $env_name

}

function install_plastome {

conda install -c bioconda getorganelle -y --name $env_name
get_organelle_config.py --add embplant_pt,embplant_mt #where are these downloaded to?
#IQTree over raxml
#IQ-TREE has a model finder to determine reconstruction parameters which RAxML didn't until recently, and still may not have it. I also like better their partition and mixture models, and less complicated command-line structure.
conda install -c bioconda iqtree -y --name $env_name
conda install -c bioconda figtree -y --name $env_name
conda install -c conda-forge mafft -y --name $env_name
}

function install_skmer {

conda install -c bioconda aliview -y --name $env_name
#seqtk
#bbmerge
#smker
#fastme
echo 'done installing packages for skmer analysis'
}

#make directories for storing results
#scripts and data directories are already made

function main {

#install_amp
#pip_install
#install_QC
#install_plastome
install_skmer
}

main


