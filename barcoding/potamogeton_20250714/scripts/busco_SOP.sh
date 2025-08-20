#!/bin/bash
#SBATCH --job-name="hibpiper"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --array=1
#SBATCH --mem=6G 

#but a few files actually doesn't exist
#well this is a very good way to turn captus into an array and run in parallel!!!

function captus_clean {

INDIR=$1
OUTDIR=$2

#set max ram =SBATCH --mem
captus clean -r $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_R*.fq.gz \
             -o $OUTDIR \
             --ram 30 \ 
             --skip_qc_stats \
             --overwrite #this will generate RESULT2 #this part is done

}

#partition=long, mem=30G
function captus_assembly {

INDIR=$1
OUTDIR=$2

#make tmp directory for megahit
mkdir $SCRATCH/tmp_${prefix}_${SLURM_ARRAY_TASK_ID}
#assembly
captus assemble -r $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_R*.fq.gz \
                -o $OUTDIR \
                --k_list 31,39,47,63,79,95,111,127 \
                --tmp_dir $SCRATCH/tmp_${prefix}_${SLURM_ARRAY_TASK_ID} \
                --overwrite \
                --ram 30  

rm -rf $SCRATCH/tmp_${prefix}_${SLURM_ARRAY_TASK_ID}  #cleaning up intermediate files

#EXPLAIN PARAMETERS:
#--k_list: adjust the kmers based on read length (150 bp)
#sample name=everything before _R*.fq.gz

#a nice command to check progress: for f in $(ls slurm-3125793_*); do echo $f; cat $f|grep 'De novo assembling with MEGAHIT:' -A 1; done
}

#extract the BUSCO genes for downstream analysis
function captus_extract {

INDIR=$1
OUTDIR=$2
REF=$3

ls $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}/01_assembly/assembly.fasta
captus extract --fasta $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}__captus-asm/01_assembly/assembly.fasta \
               --captus_assemblies_dir $OUTDIR/tmp_${prefix}_${SLURM_ARRAY_TASK_ID} \
               -n Angiosperms353 \
               --out $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID} \
               --overwrite \
               --ram 60 
               
rm -rf $OUTDIR/tmp_${prefix}_${SLURM_ARRAY_TASK_ID}

#remember to use --ram to control the amount of ram to avoid OOM error
#i don't need these. but if we want plastid and mitochondria, can be extracted as well
#--nuc_refs $REF \
#--nuc_transtable 1 \
#should i use a different translation table???
#-p SeedPlantsPTD \
#-m SeedPlantsMIT \
#don't need the per marker per sample files. just per marker all sample file for alignment is enough
}

function captus_extract_single {

INDIR=$1
OUTDIR=$2
REF=$3

captus extract --captus_assemblies_dir $INDIR \
               --out $OUTDIR/ \
               --nuc_refs $REF \
               --overwrite \
               --ram 60 \

}

function captus_extract_test {

INDIR=$1
OUTDIR=$2
REF=$3

captus extract -a $INDIR \
               --out $OUTDIR/ \
               --overwrite \
               -n Angiosperms353 \
               -c \--max_loci_files 500
}
#
function run_busco {

FASTA=$1
OUTDIR=$2
NAME=$3

busco -i $FASTA \
      --out_path $OUTDIR \
      --download_path $OUTDIR/busco_downloads \
      --augustus \
      -o $NAME \
      -c 32 \
      -m genome \
      -l eudicotyledons_odb12 \
      -f


}
#
##########################################################################################################################################
#use hybpiper instead
#run in the env hybpiper

#i think the same $BUSCO_BAITS would work here
function hyb_assemble {

INDIR=$1
OUTDIR=$2
REF=$3

echo 'hybpiper'
cd $OUTDIR

hybpiper assemble -r $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}_R*.fq.gz \
                  -t_aa $REF \
                  --prefix ${prefix}_${SLURM_ARRAY_TASK_ID} \
                  --diamond \
                  --cpu 8 \
                  --chimeric_stitched_contig_check \
                  --force_overwrite 

}

function hyb_stats {

INDIR=$1
REF=$2

cd $INDIR
ls |grep ${prefix} > namelist.txt
hybpiper stats -t_aa $REF gene namelist.txt

}

function hyb_heatmap {

INDIR=$1
REF=$2

cd $INDIR
hybpiper recovery_heatmap seq_lengths.tsv
}

#retrieve the genes
function hyb_retrieve {

INDIR=$1
REF=$2

echo 'hybpiper'
cd $INDIR
ls |grep ${prefix} > namelist.txt
hybpiper retrieve_sequences -t_aa $REF \
                            --sample_names namelist.txt \
                            --fasta_dir $contigs \
                            dna  
}

#filterings
function remove_genes_with_paralogs {

INDIR=$1
all_contigs=$2
para=$3

nb=$(cat $INDIR/${prefix}_*/*with_paralog_warnings_by_contig_depth.csv|grep True|cut -f 2 -d ','|sort|uniq|wc -l)
echo 'removing' $nb 'genes with paralog warnings'
cat $INDIR/${prefix}_*/*with_paralog_warnings_by_contig_depth.csv|grep True|cut -f 2 -d ','|sort|uniq|cut -f 3 -d ' ' > $INDIR/gene_to_remove.txt
#move the genes with paralogs to a seperate directory
for gene in  $(cat $INDIR/gene_to_remove.txt); 
do \
mv $all_contigs/${gene}.FNA $para; \
done

}

#more filtering. but still there will be about 2000 genes left. probably too many to screen?
#more filtering?
function remove_genes_with_too_few_samples {

all_contigs=$1
too_few=$2
at_least=$3

for gene in $(ls ${all_contigs}); \
do  nb=$(grep -c '>' ${all_contigs}/$gene); \
if [ $nb -lt $at_least ]; \
then  mv ${all_contigs}/$gene $too_few; \
fi;  \
done

nb=$(ls $too_few)|wc -l
echo $nb 'genes have less than' $at_least 'samples'
echo 'keep' $(ls $all_contigs|wc -l) 'genes'
}

####make alignment and phylogeny RUN IN ARRAY OF GENES!!!
function rename_contigs {

INDIR=$1
OUTDIR=$2
CSV=$3

./rename_files.py rename_contig -i $INDIR --infile gene_${SLURM_ARRAY_TASK_ID}.fasta --fcsv $CSV -o $OUTDIR
}
#
function remove_bad_samples {

INDIR=$1

./remove_taxa.py $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta
}

function RUN_mafft {

INDIR=$1
OUTDIR=$2

mafft --maxiterate 10000 $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta > $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta
}

function iqtree_per_gene {

INDIR=$1
OUTDIR=$2

mkdir $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
cp $INDIR/gene_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
#remove bad taxa (sample is bad) 'repens' and other all gap contigs
./remove_taxa.py $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta #in this case 'repens' does not exist, just remove all gap is good enough 
#remove the outliers from an alignment
./remove_highlyhetero.py $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta -t 0.4 -m majority 
iqtree -s $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}.fasta -bb 1000 -redo -safe
mv $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}/gene_${SLURM_ARRAY_TASK_ID}*treefile $OUTDIR
rm -rf $OUTDIR/gene_${SLURM_ARRAY_TASK_ID}
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

#UNIVERSAL VARIABLES=======================================================================================================================
env_name=busco
USER=zedchen 
WORKDIR=$SCRATCH/Barcoding_km/barcoding_potamogeton #
prefix=potamogeton
#SLURM_ARRAY_TASK_ID=1

####constants############
DATA=$HOME/projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/
CSV=$DATA/renamed.csv
BUSCO_BAITS=/mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/04_BUSCO/busco_downloads/lineages/eudicotyledons_odb12/ancestral_variants #from busco download
astral=$HOME/apps/manual/ASTRAL/astral.5.7.8.jar 
#==========================================================================================================================================


#subdirectories

RENAMED=$WORKDIR/results/renamed
RESULT1=$WORKDIR/results/01_clean_reads
RESULT2=$WORKDIR/results/02_captus_assembly
RESULT3=$WORKDIR/results/03_captus_extract
RESULT4=$WORKDIR/results/04_BUSCO
RESULT5=$WORKDIR/results/05_hp_assemble #hibpiper
contigs=$RESULT5/contigs_unaligned
paralog=$RESULT5/contigs_paralogs
too_few=$RESULT5/contigs_toofewsamples
RESULT6=$WORKDIR/results/06_hp_aligned
RESULT7=$WORKDIR/results/07_iqtree
RESULT8=$WORKDIR/results/08_astral

#USAGES
#=============================================
#CAPTUS: run in env captus
#captus_clean $DATA $RESULT1
#captus_assembly $RESULT1/00_adaptors_trimmed $RESULT2
#captus_extract $RESULT2 $RESULT3 $BUSCO_BAITS #memory sensiive: OOM: 15G, 30G,60G| PASS: ?G
#captus_extract_single $RESULT2 $RESULT3 $BUSCO_BAITS


#HYBPIPER ARRAY=sample size: 21-124
#hyb_assemble $DATA $RESULT5 $BUSCO_BAITS

#RUN AS SINGLE
function hybpiper_single {

#hyb_stats $RESULT5 $BUSCO_BAITS
#hyb_heatmap $RESULT5 $BUSCO_BAITS
hyb_retrieve $RESULT5 $BUSCO_BAITS
remove_genes_with_paralogs $RESULT5 $contigs $paralog
remove_genes_with_too_few_samples $contigs $too_few 100
for i in $(ls ${contigs}/*FNA); do new=${i/FNA/fasta}; mv $i $new; done
 ##just get new gene names for array

}
#hybpiper_single
#./rename_files.py rename_fasta_easy353 -i $contigs -o $contigs

#ARRAY = GENE NB: 1-709
#rename_contigs $contigs $contigs $CSV #nwo the combined samples are given different names
#remove_bad_samples $contigs
#RUN_mafft $contigs $RESULT6
#iqtree_per_gene $RESULT6 $RESULT7

#Run as single
#./count_mono.py process_treefiles -i $RESULT7 -o $RESULT8
#use astral method to make phylogeny: how many high-resolution genes per taxa do we need?
run_astral $RESULT7 $RESULT8 2 #get at least 2 genes for each mono taxa
run_astral $RESULT7 $RESULT8 3 #get at least 3 genes for each mono taxa
run_astral $RESULT7 $RESULT8 4 #get at least 4 genes for each mono taxa
run_astral $RESULT7 $RESULT8 5 #get at least 5 genes for each mono taxa
run_astral $RESULT7 $RESULT8 6 #get at least 6 genes for each mono taxa

}

main


