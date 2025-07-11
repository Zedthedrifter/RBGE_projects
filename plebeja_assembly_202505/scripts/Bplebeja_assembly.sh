#!/bin/bash
#SBATCH --job-name="Bplebeja"
#SBATCH --export=ALL
#SBATCH --partition=gpu --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G 
#Platanus takes a lot of memory

#try to run phasing step with GPU instead of CPU
#hopefully enough memory: --partition=gpu --gpus=1

#QC has been done
#remove reads mapped to chloroplast genome? I guess it doesn't matter too much
#Spade assembly

#Note, that SPAdes does not perform assembly using genomes of closely-related species. Only contigs of the same genome should be specified.
#so probably this one doesn't work
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1911-6 #a helpful paper

function quality_control {
r1=$1
r2=$2


#trim just in case
function trim {
adapters=/home/zchen/apps/conda/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa #need to specify the full path to the adapter files
trimmomatic PE $r1 $r2 ${r1/fastq/trimmed.fastq} ${r1/fastq/unpaired.fastq} ${r2/fastq/trimmed.fastq} ${r2/fastq/unpaired.fastq} \
ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
/
}
#fastqc output checking
#fastqc $r1 $r2 -o QC
trim #run once is enough
#check QC again after trimming 
fastqc ${r1/fastq/trimmed.fastq} ${r2/fastq/trimmed.fastq} -o QC
}

#===================================================================
#remove the plastmid reads
function remove_chlor {
#reference genome: NC_000932.1 Thale cress chloroplast reference genome
#update: use Vaccinium macrocarpon reference genome. a common reference for Ericaceae
ref=$1
r1=$2
r2=$3
strain=$4
indir=$5
#alignment with bowtie2
#bowtie2-build -f $ref chlo_ref
#bowtie2-align-s --wrapper basic-0 -x chlo_ref -p 20 --threads 16 --phred33 --very-sensitive --quiet --time --rg-id ${strain} --rg SM:${strain} --rg PL:'ILLUMINA' -S reads.bam -1 ${r1} -2 ${r2}
#split mapped vs unmapped reads; output the bam files into directory 'mapping'
#samtools view -b -f 4 -o ${strain}_unmapped.reads.bam reads.bam #non plasmid reads
#samtools view -b -F 4 -o ${strain}_mapped.reads.bam reads.bam #plasmid reads
#output as paired reads
#samtools fastq -1 $indir/${strain}_1.nc.fastq -2 $indir/${strain}_2.nc.fastq -0 /dev/null -s /dev/null -n ${strain}_unmapped.reads.bam
#rm ${strain}_unmapped.reads.bam
#manually gzip both outputs
gzip $indir/${strain}_1.nc.fastq
gzip $indir/${strain}_2.nc.fastq
}
#===================================================================
#map to the reference sister taxa and split the mapped/unmapped reads
#call vcf on the mapped reads
#de novo assembly for unmapped reads
function map_to_ref {
ref=$1
r1=$2
r2=$3
strain=$4
indir=$5

#alignment with bowtie2
#bowtie2-build -f $ref nc_ref
bowtie2-align-s --wrapper basic-0 -x nc_ref -p 20 --threads 16 --phred33 --very-sensitive --quiet --time --rg-id ${strain} --rg SM:${strain} --rg PL:'ILLUMINA' -S reads.bam -1 ${r1} -2 ${r2}
#split mapped vs unmapped reads; output the bam files into directory 'mapping'
samtools view -b -f 4 -o ${strain}_unmapped.reads.bam reads.bam #reads definitely for de novo assembly
samtools view -b -F 4 -o ${strain}_mapped.reads.bam reads.bam #use the mapping for variant calling?
#output as paired reads
samtools fastq -1 $indir/${strain}_1.nc.unmapped.fastq -2 $indir/${strain}_2.nc.unmapped.fastq -0 /dev/null -s /dev/null -n ${strain}_unmapped.reads.bam
gzip $indir/${strain}_1.nc.unmapped.fastq
gzip $indir/${strain}_2.nc.unmapped.fastq
rm ${strain}_unmapped.reads.bam
}
#===================================================================
#probably not a good option, but worth a try

function assembly_spade {

r1=$1
r2=$2
ref=$3
spades.py --only-assembler --trusted-contigs $ref -1 $r1 -2 $r2 -k 23,31 --careful -t 16 -m 128 -o spades_assembly
}


#===================================================================
#a de novo assembler specialized in plant genome assembly
function assembly_platanus {

r1=$1
r2=$2
prefix=$3
assembly_dir=$4

#god it's great now I've configured everything
#platanus_allee assemble -f $r1 $r2 -o $prefix -t 32 -m 196 2>assemble.log
#mv ${prefix}_contig.fa $assembly_dir

##phasing -- this module won't work on the currently available platanus version if your genome is longer than ~600Mbp. Email me for the bug-free version (lcampos@rbge.org.uk)
platanus_allee phase -o ${prefix}-phased -c $assembly_dir/${prefix}_contig.fa -IP1 $r1 $r2 -mapper minimap2 -t 24  2>phase.log #3hrs +


#platanus consensus
#platanus_allee consensus -o ${prefix}-consensus -c ${prefix}-phased_allPhaseBlock.fa -IP1 $r1 $r2 -t 72 2>consensus.log
}
#===================================================================
#scaffolding with RagTag
function Polishing_assembly {

ref=$1
query=$2
outdir=$3

ragtag.py scaffold $ref $query -o $outdir --aligner minimap2 -t 16 
#patch filling is probably introducing too much bias for similarity btw two genomes
#scaffolding is enough i think, just fill with NNNN
}
#===================================================================


#de novo assembly evaluation
function assembly_evaluation {

indir=$1
contigs=$2
outdir=$3

quast.py $indir/$contigs -o $outdir  
}
#===================================================================

function main {
name=Bplebeja
sample_dir='/home/zchen/projects/rbge/Begonia_map/Cyn1_raw/raw_data/Sample_85'
r1=$sample_dir/Sample_85_EDSW200001003-1a_H2MKTDSXY_L2_1.fq.gz
r2=$sample_dir/Sample_85_EDSW200001003-1a_H2MKTDSXY_L2_2.fq.gz
ref='/home/zchen/projects/rbge/Begonia_genomes/Reference_Assembelies/conch_genome_v4.fasta'
chlo_ref=/home/zchen/maternity_cover/Bplebeja_assembly_202505/inputs/Bplebeja_chlo.fasta
workdir='/home/zchen/maternity_cover/Bplebeja_assembly_202505'
indir=/home/zchen/maternity_cover/Bplebeja_assembly_202505/inputs
assembly_dir=/home/zchen/maternity_cover/Bplebeja_assembly_202505/assembly
#after removing reads mapped to plasmid reference:
ncr1=$indir/${name}_1.nc.trimmed.fastq.gz
ncr2=$indir/${name}_2.nc.trimmed.fastq.gz 

#QC
#quality_control $ncr1 $ncr2
#remove_chlor $chlo_ref $r1 $r2 $name $indir #done 30/05/2025
#the plasmid DNA has been assembled. we are only interested in the nuclear genome

#map to the reference genome
#map_to_ref $ref $ncr1 $ncr2 $name $indir
#busco -i $ref -o ref_BUSCO -c 32 -m genome -l eudicotyledons_odb12 -r  #check what a good genome should look like
#busco -i $assembly_dir/bp2_contig.fa -o Bp_BUSCO -c 32 -m genome -l eudicotyledons_odb12 -f
#assembly
#assembly_spade $indir/$ncr1 $indir/$ncr2 $ref
#===================================================================
assembly_platanus $ncr1 $ncr2 bp2 $assembly_dir
#Polishing_assembly $ref final_contig.fasta $assembly_dir
#===================================================================
#assembly assessment
#assembly_evaluation spades_assembly contigs.fasta spades_assembly_quast # 
#assembly_evaluation assembly bp2_contig.fa assembly_quast
}

main
