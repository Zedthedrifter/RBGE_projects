#!/bin/bash
#SBATCH --job-name="Bcassembly" 
#SBATCH --export=ALL
#SBATCH --partition=medium
#SBATCH --mem=200G 

#

function install_busco {

conda create -n Busco -c conda-forge -c bioconda busco=6.0.0 -y
#you need to manually change the !# line of busco to /home/zchen/projects/rbge/zedchen/env/Busco/bin/python3
}

function captus_assembly {

INDIR=$1

cd $INDIR/.. #got to the directory just above the result directory in use
captus clean -r $INDIR #this will generate RESULT2 #this part is done
captus assemble -r $WORKDIR/results/01_clean_reads #the general structure of captus output:RESULT3
cd - #go back to the script directory
}

function captus_extract {

FASTA=$1
REF=$2

cd $INDIR/..
captus extract -f $FASTA --dna_refs $REF
}

function quast_evaluation {

CONTIGS=$1
OUTDIR=$2

quast.py $CONTIGS -o $OUTDIR  
}

#===================================================================
#a note to yourself: if the package cannot be installed, just make another env!!!
#scaffolding with RagTag this one has to be run in a different environment called ####ragtag####!!!!! it cannot be installed with captus 
function Polishing_assembly {

INFILE=$1
OUTDIR=$2


ragtag.py scaffold $REF $INFILE -o $OUTDIR --aligner minimap2 -t 16 
#patch filling is probably introducing too much bias for similarity btw two genomes
#scaffolding is enough i think, just fill with NNNN
}

#==========================all RNA polish functions=========================================
function install_STAR {

ENV=$1

git clone https://github.com/alexdobin/STAR.git $HOME/projects/rbge/zedchen/env/$ENV/bin/STAR
cd $HOME/projects/rbge/zedchen/env/$ENV/bin/STAR/source
make STAR CXXFLAGS_SIMD=sse #some errors but seems to compile just fine
}

function RNA_polish_env_setup {

ENV=$1
#STAR: ultrafast universal RNA-seq aligner
conda create -n $ENV -y #make an env 
conda install -c bioconda samtools -y --name $ENV
conda install -c bioconda pilon -y --name $ENV
install_STAR $ENV
}

#put one fastq.gz file from each tissue type into one
#I run this in background using screen
#STEP 1
function collect_data {

INDIR=$1
OUTDIR=$2
prefix=$3

for sample in 11 45 22 65 12 #female flower, male flower, peticole, root, veg bud
do 

cat $INDIR/$sample/*_1.sanfastq.gz >> $OUTDIR/${prefix}_RNA_R1.fastq.gz
cat $INDIR/$sample/*_2.sanfastq.gz >> $OUTDIR/${prefix}_RNA_R2.fastq.gz

done
}

#STEP 2
#Align RNA-seq Reads to the Genome
function STAR_index {

OUTDIR=$1
genome=$2
prefix=$3

cd $OUTDIR
echo 'Running STAR'
rm $OUTDIR/genomeDir -rf
mkdir $OUTDIR/genomeDir
chmod 777 $OUTDIR/genomeDir

#make a genome index
$STAR --runThreadN 16 \
     --genomeDir $OUTDIR/genomeDir \
     --runMode genomeGenerate \
     --genomeFastaFiles $genome \
     --genomeChrBinNbits  18 \
     --genomeSAindexNbases 13 \ 
     --limitGenomeGenerateRAM 60000000000 
     #60G RAM limit
 

}

#STEP 3
#STAR MAPPING
function STAR_map {

OUTDIR=$1
genome=$2
prefix=$3

cd $OUTDIR
echo 'Running STAR'

#mapping
$STAR --runThreadN 16 \
     --genomeDir $OUTDIR/genomeDir \
     --readFilesIn $OUTDIR/${prefix}_RNA_R1.fastq.gz $OUTDIR/${prefix}_RNA_R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix $OUTDIR/Bp 
     #specify the output directory as well
}

#STEP 4.0
function samtool_index {

OUTDIR=$1
genome=$2
prefix=$3

samtools view -Sbh -F 4 $OUTDIR/${prefix}Aligned.sortedByCoord.out.bam|samtools sort -o $OUTDIR/${prefix}_mapped.sorted.bam
samtools index $OUTDIR/${prefix}_mapped.sorted.bam
}

#STEP 4.1
#improve assembly with Pilon: output the corrected fasta file based on bam mapping
#the source code of the executable has to be changed, otherwise default max RAM=1G is way too small
#100G: OOM; 200 G: worked 
function Pilon_polish {

OUTDIR=$1
genome=$2
prefix=$3

pilon --genome $genome \
      --frags $OUTDIR/${prefix}_mapped.sorted.bam \
      --output ${prefix}_pilon \
      --outdir $OUTDIR \
      --vcf \
      --changes

}

#STEP 5.0
#INSTALL BRAKER in a new ENV
function Braker_installation {

conda create -n braker -y
#install dependencies manually
conda install -c anaconda perl -y -n braker
conda install -c anaconda biopython -y -n braker
conda install -c bioconda perl-app-cpanminus -y -n braker
cpanm --local-lib=~/perl5 YAML
conda install -c bioconda perl-file-spec -y -n braker
conda install -c bioconda perl-hash-merge -y -n braker
conda install -c bioconda perl-list-util -y -n braker
conda install -c bioconda perl-module-load-conditional -y -n braker
conda install -c bioconda perl-posix -y -n braker
conda install -c bioconda perl-file-homedir -y -n braker
conda install -c bioconda perl-parallel-forkmanager -y -n braker
conda install -c bioconda perl-scalar-util-numeric -y -n braker
conda install -c bioconda perl-yaml -y -n braker
conda install -c bioconda perl-class-data-inheritable -y -n braker
conda install -c bioconda perl-exception-class -y -n braker
conda install -c bioconda perl-test-pod -y -n braker
conda install -c bioconda perl-mce -y -n braker
conda install -c bioconda perl-threaded -y -n braker
conda install -c bioconda perl-list-util -y -n braker
conda install -c bioconda perl-math-utils -y -n braker
conda install -c bioconda cdbtools -y -n braker
conda install -c eumetsat perl-yaml-xs -y -n braker
conda install -c bioconda perl-data-dumper -y -n braker
conda install -c bioconda perl-Logger-Simple -y -n braker
conda install -c bioconda samtools -y -n braker
conda install -c bioconda bamtools -y -n braker
conda install -c bioconda augustus -y -n braker
#other required tools for braker
conda install -c bioconda diamond -y -n braker
conda install -c bioconda bedtools -y -n braker
conda install -c bioconda gffread -y -n braker
#finally, installing braker
conda install -c bioconda braker  -y -n braker
}
#mandatory tools for braker not available on conda 
function stringtie2_installation {

git clone https://github.com/gpertea/stringtie $HOME/projects/rbge/zedchen/manual/stringtie
cd $HOME/projects/rbge/zedchen/manual/stringtie
make release
}

function genemark_installation {

git clone https://github.com/gatech-genemark/GeneMark-ETP  $HOME/projects/rbge/zedchen/manual/genemark
cd $HOME/projects/rbge/zedchen/manual/genemark
./check_install.pl
}

function braker_config {
#yeah this is a good idea...keeping configuration within the script
#so nothing happening within a certain env will leak to the global config
export PATH=/$HOME//projects/rbge/zedchen/manual/stringtie:$PATH
export AUGUSTUS_CONFIG_PATH=/home/zchen/projects/rbge/zedchen/pkg/augustus-3.1-0/scripts
export GENEMARK_PATH=$HOME/projects/rbge/zedchen/manual/genemark/bin/gmes #in a subdirectory called gmes
export PERL5LIB=$HOME/projects/rbge/zedchen/manual/perl5/lib/perl5:$PERL5LIB #configure to god knows where and why
export PERL5LIB=$HOME/projects/rbge/zedchen/env/braker/lib/site_perl/5.26.2/:$PERL5LIB #configure to env
}
#STEP 5.1
#Gene prediction using braker for the polished genome #also a way to test the quality of assembly
function Braker_polish {

OUTDIR=$1
genome=$2
prefix=$3
species=$4


braker.pl --genome=${genome} \
          --bam=$OUTDIR/${prefix}_mapped.sorted.bam \
          --BAMTOOLS_PATH=/home/zchen/projects/rbge/zedchen/env/braker/bin/ \
          --species=$species \
          --gff3 \
          --overwrite \
          --workingdir=$OUTDIR 
          
}


#MASTER function
function RNA_polish_assembly {

GENDIR=$1
genome=$2
prefix=$3

#RNA_polish_env_setup RNA_polish

#collect_data $HOME/projects/rbge/Begonia_RNASeq/raw-reads/Con_Ple_Raw_KE $GENDIR $prefix

#STAR_index $GENDIR $genome $prefix #done, didn't take too long
#STAR_map $GENDIR $genome $prefix

#samtool_index $GENDIR $genome $prefix

#mem=200G
#Pilon_polish $GENDIR $genome $prefix #need to change max Ram in the executable. stupid default = 1G

#ALL THE PREP WORK BEFORE RUNNING BRAKER (in a different env)
#Braker_installation #only need to do once. now we can use env: braker
#stringtie2_installation
#genemark_installation

#Config in ENV for braker
braker_config

#god i can't belive how much preparation it took to run braker...
#mem=200G
Braker_polish $GENDIR $GENDIR/${prefix}_pilon.fasta $prefix Begonia_plebeja #run in ENV braker (installation incompatibility with other packages)
}

#==============================END OF RNA POLISHING FUNCTIONS=====================================


#===================================================================

function main {

#conda create -n captus -c bioconda -c conda-forge captus iqtree -y# T

#UNIVERSAL VARIABLES
env_name=captus 
USER=zedchen 
WORKDIR=${SCRATCH}/plebeja_assembly_202505
prefix=Bp

#CONSTANTS
REF=$HOME/projects/rbge/Begonia_genomes/Reference_Assembelies/conch_genome_v4.fasta
pip=$HOME/projects/rbge/$USER/env/easy353/bin/pip
STAR=$HOME/projects/rbge/zedchen/env/RNA_polish/bin/STAR/bin/Linux_x86_64/STAR
#subdirectories
DATA=$WORKDIR/inputs
#RNA=
RESULT1=$WORKDIR/results/00_qctrimm 
RESULT2=$WORKDIR/results/01_clean_reads #created by captus
RESULT3=$WORKDIR/results/02_assemblies #created by captus
RESULT4=$WORKDIR/results/03_extractions
RESULT5=$WORKDIR/results/04_quast_evaluation
RESULT6=$WORKDIR/results/05_ragtagscaffolding
RESULT7=$WORKDIR/results/06_BUSCO
RESULT8=$WORKDIR/results/07_RNA_polish
#$WORKDIR/results/
#mkdir $RESULT3

#USAGES
#=============================================
#captus_assembly $RESULT1 #took two days, must be submitted as a long job

#scaffolding #oh this works very well!!
#mkdir $RESULT6/02_assemblies
#Polishing_assembly $RESULT3/Bplebeja.trimmed__captus-asm/01_assembly/assembly.fasta $RESULT6/02_assemblies #run on env ragtag!

#the size of the assembly is much larger than the ref, use extraction to select for homologous genome?
#might as well write the commands directly. you have to specify the fasta file when using non-captus assembly
cd $WORKDIR/results
#captus extract -a $RESULT4 -f $RESULT6/02_assemblies/ragtag.scaffold.fasta --dna_refs $REF #reference based assembly
#extracting after scaffolding takes days,just do the extraction before scaffolding
#captus extract -a $RESULT4 -f $RESULT3/Bplebeja.trimmed__captus-asm/01_assembly/assembly.fasta --dna_refs $REF #doesn't work, OOM, just extract from ragtag results
#Polishing_assembly $RESULT3/Bplebeja.trimmed__captus-asm/01_assembly/assembly.fasta $RESULT6/02_assemblies #run on env ragtag!
cd -

#Polish with RNA-seq data:
#input variables: OUTDIR genome prefix
RNA_polish_assembly $RESULT8 $RESULT6/ragtag.scaffold.only.fasta Bp #polish the scaffolded assembly with transcriptom


#quality check: completeness of gene space: BUSCO 
#run in conda env busco
#install_busco
cd $RESULT7 #run busco in the same directory
#check the completeness of the reference genome
#busco -i $REF -o REF -c 32 -m genome -l eudicotyledons_odb12 -r  #check what a good genome should look like
#directly from scaffolded assembly
#busco -i $RESULT6/02_assemblies/ragtag.scaffold.fasta -o ragtag -c 32 -m genome -l eudicotyledons_odb12 -f
#the extraction results
#busco -i $RESULT4/?????? -o $RESULT7/extract -c 32 -m genome -l eudicotyledons_odb12 -f
cd - #come back
}

main


