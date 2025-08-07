#DNA barcoding project, salix

srun --partition=short --cpus-per-task=8 --mem=16G --pty bash

create env salix

conda activate salix

./setup_env_salix #this can be resued later for potamogeton

#rename, QC and cleaning up reads: arrat
sbatch cleaning_QC.sh

#plastome phylgeny
download reference genome: Salix blinii
https://www.ncbi.nlm.nih.gov/nuccore/PP833162.1

sbatch plastome_phylogeny.sh (mafft, iqtree)

#skmer

#Angiosperm353

single: setup_env
single: build_353_db
array: easy353_assembly
single:python $easy353path/script/multi_sample_gene_merge.py
single:python $easy353path/script/correct_seq_ori.py
single:./rename_files.py rename_fasta_easy353  #rename gene assemlby so that we can prepare for mafft using array
gene_array:easy353_mafft
gene_array: rename contigs
gene_array:easy353_trimal #filter the alignments
#get a summary of alignment: a python script
./rename_files.py gene_alignment_summary -i $RESULT9 -o $RESULT9
#filter sample and genes and concatenate the selected aligned contigs
./rename_files.py easy353_salix_filter -i $RESULT9 -o $RESULT10
filtered_array: easy353_trimal $RESULT10 $RESULT10 #another trimming with trimal
#Astral method


#Concatenated
filtered_array: easy353_phylogeny $RESULT10 $RESULT11#build phylogeny with IQtree
#BUSCO

#Genome-wide SNPs 