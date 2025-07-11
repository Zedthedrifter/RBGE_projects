#!/usr/bin/env bash
#SBATCH --partition=long 	# 
#SBATCH --time=0-72:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12

## compressing raw files

#cd /data/camposl/dregei_scaffolding/

#gzip ./dregei_forward_paired.fq
#gzip ./dregei_reverse_paired.fq

##platanus on dregei

cd /mnt/shared/scratch/lcampos/dregei_assembly/

##contig assembly

/mnt/shared/scratch/lcampos/software/Platanus_allee_v2.2.2_Linux_x86_64/./platanus_allee assemble \
-f dregei_forward_paired.fq.gz dregei_reverse_paired.fq.gz \
-o dregei-platanus-contigs -t 48 -m 128 2>assembly.log

##phasing -- this module won't work on the currently available platanus version if your genome is longer than ~600Mbp. Email me for the bug-free version (lcampos@rbge.org.uk)

/mnt/shared/scratch/lcampos/software/Platanus_allee_v2.2.2_Linux_x86_64/./platanus_allee phase \
-o dregei-platanus-phased \
-c ./dregei-platanus-contigs_contig.fa \
-IP1 dregei_forward_paired.fq.gz dregei_reverse_paired.fq.gz \
-mapper /mnt/shared/scratch/lcampos/software/Platanus_allee_v2.2.2_Linux_x86_64/minimap2-2.0-r191/./minimap2 \
-t 32  2>phase.log

#platanus consensus

/mnt/shared/scratch/lcampos/software/Platanus_allee_v2.2.2_Linux_x86_64/./platanus_allee consensus \
-o dregei-platanus-consensus \
-c ./dregei-platanus-phased_consensusInput.fa \
-IP1 luxurians_forward_paired.fq.gz luxurians_reverse_paired.fq.gz \
-t 72 2>consensus.log

exit
