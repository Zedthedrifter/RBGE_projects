#!/bin/bash
#SBATCH --job-name="vcffasta"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --array=1
#SBATCH --mem=2G

#run MuMmer
function chromosome_homolog {

reference=$1
query=$2
OUTDIR=$3

mkdir $OUTDIR/chromosome${SLURM_ARRAY_TASK_ID}
cd $OUTDIR/chromosome${SLURM_ARRAY_TASK_ID}
#dnadiff is in captus env
dnadiff $reference/Sh_chromosome_${SLURM_ARRAY_TASK_ID}.fasta \
        $query/chromosome_${SLURM_ARRAY_TASK_ID}.fasta

}

#STEP1: Indexing
#ENV: snps
function chrom_index {

REFDIR=$1
FASTA=$2
IDX=$3

cd $REFDIR
bowtie2-build -f $FASTA $IDX
}

#STEP2: MAPPING
function chrom_map {

INDIR=$1
IDX=$2
OUTDIR=$3

echo 'RUNNING BOWTIE2 MAPPING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
bowtie2-align-s --wrapper basic-0 \
                -x $IDX \
                -1 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R1.fq.gz \
                -2 ${INDIR}/${prefix}_${SLURM_ARRAY_TASK_ID}_R2.fq.gz \
                -p 20 \
                -N 1 \
                -L 20 \
                --threads 8 \
                --phred33 \
                --very-sensitive \
                --no-discordant \
                --no-mixed \
                --no-unal \
                --time \
                --rg-id ${prefix}_${SLURM_ARRAY_TASK_ID} \
                --rg SM:${prefix}_${SLURM_ARRAY_TASK_ID} \
                --rg PL:'ILLUMINA' |\
                samtools view -Sbh -F 4 -f 3 -q 30 -@ 8|\
                samtools sort -@ 8 -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

echo 'RUNNING SAMTOOLS COVERAGE ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools coverage $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam > $SORTED/${prefix}_${SLURM_ARRAY_TASK_ID}_depth.txt
}
#a patch
function index_bam {

OUTDIR=$1

echo 'RUNNING SAMTOOLS INDEXING ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools index $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam

echo 'RUNNING SAMTOOLS COVERAGE ON' ${prefix} ${SLURM_ARRAY_TASK_ID}
samtools coverage $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam > $SORTED/${prefix}_${SLURM_ARRAY_TASK_ID}_depth.txt

}

#need to do this for sample1-16 once. i can do this in a simple for loop though; array is so much faster
function fix_RG_bam {

OUTDIR=$1

RG='@RG\tID:'${prefix}_${SLURM_ARRAY_TASK_ID}'\tSM:'${prefix}_${SLURM_ARRAY_TASK_ID}'\tPL:ILLUMINA'

samtools addreplacerg -r $RG \
                      -o $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.fixed.bam \
                      $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.sorted.bam
}

#STEP3: SNP calling : 
function snp_calling {

INDIR=$1
OUTDIR=$2
REF=$3
NAME=$4

echo 'calling SNPs'
ls $INDIR/${prefix}*sorted.bam > bamlist.txt #need to change back to sorted.bam

bcftools mpileup -Ou -f $REF --bam-list bamlist.txt --threads 80| \
bcftools call -Ou -mv | \
bcftools filter -s LowQual -e 'QUAL<20 || DP<10' > $OUTDIR/${NAME}.flt1.vcf #combined depth across samples>100 or quality >20
#this is a soft filter, which does not remove any snps. you can remove them later in the df easily

#rm bamlist.txt

}

#STEP4: SNP filtering
function filter_vcf {

INDIR=$1
OUTDIR=$2
NAME=$3

# Prep: Filter SNPs using vcftools

vcftools --vcf $INDIR/${NAME}.flt1.vcf \
         --out $OUTDIR/${NAME}.flt2 \
         --recode --recode-INFO-all \
         --minQ 30 \
         --max-missing 0.95 \
         --minDP 100 \
         --maxDP 100000 \
         --min-alleles 2 \
         --hwe 0.05

#minDP=100, 2x per sample (47). very fair here
#maxDP=47x250=11750. honestly this filter is unnecessary
#I removed --mac 5, which requires allele count to be >=5. Allele count is simply the number of times that allele appears over all individuals at that site. this varies for sample size etc. 
#--max-alleles 2: no max allele
#--remove-indels: keep indels? you can remove it from the csv later
}

function rename_contigs {

INDIR=$1
OUTDIR=$2
CSV=$3

./rename_files.py rename_contig -i $INDIR --infile gene_${SLURM_ARRAY_TASK_ID}.fasta --fcsv $CSV -o $OUTDIR
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

#EXECUTION=================================================================================================================================

function main {

#UNIVERSAL VARIABLES=======================================================================================================================
env_name=captus #dnadiff is in captus.bowtie2: easy353; snps: after mapping
USER=zedchen 
WORKDIR=$SCRATCH/Barcoding_km/SNP_salix 
prefix=salix
#SLURM_ARRAY_TASK_ID=1

####constants############
DATA=$HOME/projects/rbge/pholling/barcoding/salix_combined/results/00_reads/
CSV=/mnt/shared/projects/rbge/pholling/barcoding/salix_combined/results/00_reads/renamed.csv
bcftools=$HOME/apps/manual/bcftools-1.22/bcftools #somehow it's broken and i had to install bcftools from source
astral=$HOME/apps/manual/ASTRAL/astral.5.7.8.jar 
REF=$WORKDIR/refs
Scin=$REF/S_cinerea.chromosome.fasta
Sher=$REF/S_herbacea.chromosome.fasta
Scap=$REF/S_caprea.chromosome.fasta
Sret=$REF/S_reticulata.chromosome.fasta
#==========================================================================================================================================


#subdirectories
REFR1=$WORKDIR/results/ref_01_mummer #REF RESULT
RESULT1=$WORKDIR/results/01_sorted_bam #Skip the sam stage--> pipe to bam directly
RESULT2=$WORKDIR/results/02_VCFs
RESULT3=$WORKDIR/results/03_SNP_fasta
RESULT4=$WORKDIR/results/04_IQTREE
RESULT5=$WORKDIR/results/05_ASTRAL
RESULT=$WORKDIR/results/0

#USAGES
#=============================================
function setup_dir {

mkdir $WORKDIR
mkdir $WORKDIR/results
mkdir $WORKDIR/refs
mkdir $REFR1
mkdir $RESULT1
mkdir $RESULT1/CIN
mkdir $RESULT2
mkdir $RESULT2/CIN
mkdir $RESULT3
mkdir $RESULT3/CIN160
mkdir $RESULT3/CIN150
mkdir $RESULT3/CIN140_160
mkdir $RESULT3/CIN40_50
mkdir $RESULT3/CIN160_renamed
mkdir $RESULT3/CIN150_renamed
mkdir $RESULT3/CIN140_160_renamed
mkdir $RESULT3/CIN40_50_renamed
mkdir $RESULT4
mkdir $RESULT4/CIN160
mkdir $RESULT4/CIN150
mkdir $RESULT4/CIN140_160
mkdir $RESULT5
mkdir $RESULT5/CIN160
mkdir $RESULT5/CIN150 
mkdir $RESULT5/CIN140_160 
}
setup_dir #WILL NOT OVERWRITE DIR, JUST LEAVE IT

#REFERENCE CURATION
#./process_ref.py $REF $REF
#chromosome_homolog: 10G mem
#chromosome_homolog $WORKDIR/refs  $WORKDIR/refs $REFR1 #

#STEP_1: INDEXING: SINGLE
function index_ref {

chrom_index $REF $Scin CIN
chrom_index $REF $Sher HER
chrom_index $REF $Scap CAP
chrom_index $REF $Sret RET
}
#index_ref

#STEP_2: MAPPING: ARRAY-SAMPLE 1-47
#mem=9G is enough. none exceeded 8G when processing the big salix 47 files (must be less for other salix/potamogeton)
#chrom_map $DATA $REF/CIN $RESULT1/CIN 
#index_bam #just to fix 1-13. there was a bug and they could not be indexed
#fix_RG_bam $RESULT1/CIN 

#STEP_3: SNP CALLING: SINGLE
#snp_calling $RESULT1/CIN $RESULT2/CIN $Scin Scin
#filter_vcf $RESULT2/CIN $RESULT2/CIN Scin

#STEP_4: convert VCF to fasta: SINGLE
#probably can change the names here
#collect the chromosome info into  csv
function extract_info {

INDIR=$1

rm $INDIR/clip_positions.csv -f
for fa in $(ls $INDIR/*fasta); do a=$(head $fa -n1|cut -d ' ' -f 2) ; b=$(echo $fa|rev|cut -d '/' -f 1|rev); echo $b , $a >> $INDIR/clip_positions.csv; done

}

#./vcf_csv_fa.py --vcf $RESULT2/CIN/Scin.flt2.recode.vcf --workdir $RESULT3/CIN160 --minlen 160 --maxlen 200 #
#extract_info $RESULT3/CIN160
#./vcf_csv_fa.py --vcf $RESULT2/CIN/Scin.flt2.recode.vcf --workdir $RESULT3/CIN150 --minlen 150 --maxlen 200 #634, probably the best
#extract_info $RESULT3/CIN150
#./vcf_csv_fa.py --vcf $RESULT2/CIN/Scin.flt2.recode.vcf --workdir $RESULT3/CIN140_160 --minlen 140 --maxlen 160 #way too many clips selected if you set a wide range
#extract_info $RESULT3/CIN140_160
#./vcf_csv_fa.py --vcf $RESULT2/CIN/Scin.flt2.recode.vcf --workdir $RESULT3/CIN40_50 --minlen 40 --maxlen 50 #try SNPs on a more conserved region
#extract_info $RESULT3/CIN40_50

#STEP_5: RENAME
#./rename_files.py rename_fasta_easy353 -i $RESULT3/CIN160 -o $RESULT3/CIN160_renamed
#./rename_files.py rename_fasta_easy353 -i $RESULT3/CIN150 -o $RESULT3/CIN150_renamed
#./rename_files.py rename_fasta_easy353 -i $RESULT3/CIN140_160 -o $RESULT3/CIN140_160_renamed
#./rename_files.py rename_fasta_easy353 -i $RESULT3/CIN40_50 -o $RESULT3/CIN40_50_renamed

#STEP_6:IQTREE PER GENE:gene array: env=easy353 MEM=1G
#rename_contigs $RESULT3/CIN160_renamed $RESULT3/CIN160_renamed $CSV #1-308


#iqtree_per_gene $RESULT3/CIN160_renamed $RESULT4/CIN160 #1-308
#rename_contigs $RESULT3/CIN150_renamed $RESULT3/CIN150_renamed $CSV #1-825
#iqtree_per_gene $RESULT3/CIN150_renamed $RESULT4/CIN150 #1-825 
#rename_contigs $RESULT3/CIN140_160_renamed $RESULT3/CIN140_160_renamed $CSV #1-978
#iqtree_per_gene $RESULT3/CIN140_160_renamed $RESULT4/CIN140_160 #1-978

#STEP_7: ASTRAL TREE: Run as single
#./count_mono.py process_treefiles -i $RESULT4/CIN150 -o $RESULT5/CIN150
#./count_mono.py process_treefiles -i $RESULT4/CIN160 -o $RESULT5/CIN160
#./count_mono.py process_treefiles -i $RESULT4/CIN140_160 -o $RESULT5/CIN140_160
#use astral method to make phylogeny: how many high-resolution genes per taxa do we need?
function RUN_ASTRAL {

INDIR=$1
OUTDIR=$2

run_astral $INDIR $OUTDIR 2 #get at least 2 genes for each mono taxa
run_astral $INDIR $OUTDIR 3 #get at least 3 genes for each mono taxa
run_astral $INDIR $OUTDIR 4 #get at least 4 genes for each mono taxa
}

#RUN_ASTRAL $RESULT4/CIN150 $RESULT5/CIN150
#RUN_ASTRAL $RESULT4/CIN160 $RESULT5/CIN160
RUN_ASTRAL $RESULT4/CIN140_160 $RESULT5/CIN140_160

#should i look for species specific genotypes?
#or just turn vcf to fasta --> chop --> filter --> phylogeny --> astral
 
}

main