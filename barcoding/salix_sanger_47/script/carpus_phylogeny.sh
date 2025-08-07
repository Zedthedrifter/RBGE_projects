#!/bin/bash
#SBATCH --job-name="amplicon_seq_GPU"
#SBATCH --export=ALL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1 #ask for one. without this line I'll get this error message: [error] CUDA device requested but no devices found.
#SBATCH --mem=32G #CHANGE TO THE AMOUNT OF MEM YOU NEED

salix_fq=/mnt/shared/projects/rbge/willow/new_seq_data/