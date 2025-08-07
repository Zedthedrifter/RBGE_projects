#DNA barcoding project, salix old 47 samples from sanger sequencing


conda activate salix #use this one for cleaning up data

#rename, QC and cleaning up reads: arrat
sbatch cleaning_QC.sh
  rename, use the renamed csv to map back the species

#plastome phylgeny

#skmer

#Angiosperm353

single: setup_env
single: build_353_db
array:fastp
      cp the fastp results to projects for safe keeping
array: easy353_assembly
single:python $easy353path/script/multi_sample_gene_merge.py
single:python $easy353path/script/correct_seq_ori.py
single:./rename_files.py rename_fasta_easy353  #rename gene assemlby so that we can prepare for mafft using array
#after merging: 353 genes
gene_array:easy353_mafft
gene_array: rename contigs
gene_array:easy353_trimal #filter the alignments
#get a summary of alignment: a python script
./rename_files.py gene_alignment_summary -i $RESULT9 -o $RESULT9
#filter sample and genes and concatenate the selected aligned contigs
./rename_files.py easy353_salix_filter -i $RESULT9 -o $RESULT10
filtered_array: easy353_trimal $RESULT10 $RESULT10 #another trimming with trimal
filtered_array: easy353_phylogeny $RESULT10 $RESULT11#build phylogeny with IQtree
#astral method/concatenation of the best genes
gene_array: iqtree_per_gene
single:./count_mono.py#concatenation method
#BUSCO

#Genome-wide SNPs 