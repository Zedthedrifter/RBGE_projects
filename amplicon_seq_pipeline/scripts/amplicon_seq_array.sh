#!/bin/bash
#SBATCH --job-name="ampsort"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --array=1-48 #
#SBATCH --mem=1M #change it to the amount of memory you need

SLURM_CPUS_PER_TASK=8

#parallel
function samtools_fastq {
echo 'running samtools for ${work_dir}'
INDIR=$1 
OUTDIR=$2
FILE_PREFIX=$3

#SLURM_ARRAY_TASK_ID=$(printf "%02d" "$((10#$SLURM_ARRAY_TASK_ID))")

#get all the rests, regardless of pairing etc.
samtools fastq -@ $SLURM_CPUS_PER_TASK $INDIR/${FILE_PREFIX}${SLURM_ARRAY_TASK_ID}.bam > $OUTDIR/samtools_${SLURM_ARRAY_TASK_ID}.fastq

echo 'samtools_fastq Job completed'
}


#probably not necessary most time but this time we need to merge the files
function merge_fastqgz {

INDIR=$1
OUTDIR=$2
prefix=$3
#SLURM_ARRAY_TASK_ID=$(printf "%02d" "$((10#$SLURM_ARRAY_TASK_ID))")

zcat $INDIR/AYZ217_pass_barcode${SLURM_ARRAY_TASK_ID}* > $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.fastq

}


#parallel
function NanoPlot_QC {

INDIR=$1 
OUTDIR=$2
minlen=$3
maxlen=$4
prefix=$5

#SLURM_ARRAY_TASK_ID=$(printf "%02d" "$((10#$SLURM_ARRAY_TASK_ID))")

echo 'running NanoPlot QC' $OUTDIR $minlen $maxlen

NanoPlot -t 2 --fastq $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.fastq \
  --outdir $OUTDIR/NanoPlot_bc${SLURM_ARRAY_TASK_ID}  \
  --minlength $minlen --maxlength $maxlen --plots dot --legacy hex

echo 'NanoPlot QC Job completed'

}

#parallel
function Porechop_trim {

INDIR=$1 
OUTDIR=$2
prefix=$3

#SLURM_ARRAY_TASK_ID=$(printf "%02d" "$((10#$SLURM_ARRAY_TASK_ID))")

echo 'running Porechop to trim barcode and adaptor'

porechop -t 4 --extra_end_trim 0 -i $INDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.fastq -o $OUTDIR/porechop.bc${SLURM_ARRAY_TASK_ID}.fastq

echo 'Porechop_trim Job completed'

}

#parallel
function Chopper_QT { 

INDIR=$1 
OUTDIR=$2

#SLURM_ARRAY_TASK_ID=$(printf "%02d" "$((10#$SLURM_ARRAY_TASK_ID))")

echo 'running chopper for quality trimming'

chopper -q 10 -l 500 -i $INDIR/porechop.bc${SLURM_ARRAY_TASK_ID}.fastq > $OUTDIR/chopper.bc${SLURM_ARRAY_TASK_ID}.fastq

echo 'Chopper Job completed'
}

#parallel
function ampliconsorter {

INDIR=$1 
OUTDIR=$2

#SLURM_ARRAY_TASK_ID=$(printf "%02d" "$((10#$SLURM_ARRAY_TASK_ID))")


echo 'running ampliconsorter to create consensus'
ls $INDIR/chopper.bc${SLURM_ARRAY_TASK_ID}.fastq
$PYTHON3 $AMPLICON_SORTER_PATH/amplicon_sorter.py \
        -i $INDIR/chopper.bc${SLURM_ARRAY_TASK_ID}.fastq \
        -o $OUTDIR/ampsorter.bc${SLURM_ARRAY_TASK_ID}.def \
        -np $SLURM_CPUS_PER_TASK \
        -ar -maxr 570000

echo 'ampliconsorter Job completed'
}

function setup_workdir {

mkdir $WORKDIR -p
mkdir $RESULT0 -p
mkdir $RESULT1 -p
mkdir $RESULT2 -p
mkdir $RESULT3 -p
mkdir $RESULT4 -p
mkdir $RESULT5 -p
mkdir $RESULT6 -p
}

# Check command-line argument and call the function
function main {

USER=$1 #
env_name=$2 #same as your previous input
WORKDIR=$HOME/projects/rbge/zedchen/Flongle/0903 #on scratch $SCRATCH/FLongle_0903 #
FASTQ=$HOME/projects/rbge/jal/Flongle/flongel6_data
prefix=barcode
#SLURM_ARRAY_TASK_ID=1
#---------------------------------------
SLURM_ARRAY_TASK_ID=$(printf "%02d" "$((10#$SLURM_ARRAY_TASK_ID))")
RESULT0=$WORKDIR/00-mergefq
RESULT1=$WORKDIR/01-samtools-fastq
RESULT2=$WORKDIR/02-nanoplot
RESULT3=$WORKDIR/03-porechop
RESULT4=$WORKDIR/04-chopper
RESULT5=$WORKDIR/05-ampliconsorter
RESULT6=$WORKDIR/06-consensus

#set up directory
setup_workdir

#merge_fastqgz $FASTQ $RESULT0 $prefix
#samtools 
#samtools_fastq $RESULT0 $RESULT1 ${FILE_PREFIX}
#NanoPlot_QC
min_len=$4
max_len=$5
#NanoPlot_QC $FASTQ $RESULT2 $min_len $max_len $prefix
#Porechop
#Porechop_trim $RESULT0 $RESULT3 $prefix
#Chopper
#Chopper_QT $RESULT3 $RESULT4
#ampliconsorter
#ampliconsorter $RESULT4 $RESULT5
#collect consensus: run just once
mv $RESULT5/ampsorter.bc${SLURM_ARRAY_TASK_ID}.def/consensusfile.fasta $RESULT6/consensus_${SLURM_ARRAY_TASK_ID}.fasta
}


#add PATH to dorado
export PATH=$PATH:$HOME/scratch/apps/dorado-0.4.1-linux-x64/bin
export PATH=$PATH:$HOME/scratch/apps/dorado-0.4.1-linux-x64
export LD_LIBRARY_PATH=$HOME/scratch/apps/dorado-0.4.1-linux-x64/lib:$LD_LIBRARY_PATH
MODEL=$HOME/scratch/apps/dorado-0.4.1-linux-x64/model/dna_r10.4.1_e8.2_400bps_hac\@v4.2.0/
#CONSTANTS
PYTHON3=$HOME/apps/env/minion/bin/python3
AMPLICON_SORTER_PATH=$HOME/apps/manual/amplicon_sorter
#---------------------------------------------------------------------------------------------------------------------------------
#USER INPUTS (CAN ALSO MAKE IT INTERACTIVE, LIKE INPUT FROM TERMINAL INSTEAD OF IN THE SCRIPT, IF YOU PREFER)
USER=zedchen
ENV=minion
FILE_PREFIX=SQK-NBD114-96_barcode
min_len=200
max_len=10000


main $USER $ENV $min_len $max_len #CHANGE BOTH TO YOUR ACTUAL USER NAME, ENV NAME, AND WORK DIRECTORY ON SCRATCH





