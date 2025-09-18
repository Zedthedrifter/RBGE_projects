#!/bin/bash

#this is a one time script. i used it to merge the reads from the same samples in batch 3. in the early batch, we have some xxa and xxb files, so they are merged into xx as well
#it takes some careful coordination to keep the numbers right... or at least now i believe the numbers are still accurate
#if there is any weird assembly results...just redo the merging from the source. no other way 

function main {

INDIR1=/mnt/shared/projects/rbge/pholling/barcoding/salix_combined/results/batch3_reads/
INDIR2=$SCRATCH/Barcoding_km/salix_20250910/01_fastp/batch3
OUTDIR=tmp

mkdir $OUTDIR

function merge2 {

N1=$1
N2=$2

cat $INDIR1/salix_${N1}_R1.fq.gz $INDIR2/salix_${N2}_R1.fq.gz > $OUTDIR/salix_${N2}_R1.fq.gz
cat $INDIR1/salix_${N1}_R2.fq.gz $INDIR2/salix_${N2}_R2.fq.gz > $OUTDIR/salix_${N2}_R2.fq.gz

}

function merge3 {

#TWO FILES FROM INDIR1 ONE FROM INDIR2

N1=$1
N2=$2
N3=$3

cat $INDIR1/salix_${N1}_R1.fq.gz $INDIR1/salix_${N2}_R1.fq.gz $INDIR2/salix_${N3}_R1.fq.gz > $OUTDIR/salix_${N3}_R1.fq.gz
cat $INDIR1/salix_${N1}_R2.fq.gz $INDIR1/salix_${N2}_R2.fq.gz $INDIR2/salix_${N3}_R2.fq.gz > $OUTDIR/salix_${N3}_R2.fq.gz
}

#merge2 76 76
#merge2 77 77
#merge2 78 78
#merge2 79 79
#merge2 80 80
#merge2 83 82
#merge2 84 83
#merge2 87 85
#merge2 88 86

merge3 81 82 81
merge3 85 86 84
}

main