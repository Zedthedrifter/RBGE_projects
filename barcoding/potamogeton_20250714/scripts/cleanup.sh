#!/bin/bash
#SBATCH --job-name="Potaclean"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=8
#SBATCH --mem=1G

#curation of the potamogeton data based on James' comment 2025/8/13

function remove_combine {

#remove undetermined reads (just reads from all samples unassigned)
rm projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_124*
#merge berchtoldi5 and berchtoldi6
cat projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_110_R1.fq.gz \
    projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_111_R1.fq.gz \
    > projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_125_R1.fq.gz

cat projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_110_R2.fq.gz \
    projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_111_R2.fq.gz \
    > projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_125_R2.fq.gz

#merge filiformis4 and filimormis5
cat projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_112_R1.fq.gz \
    projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_113_R1.fq.gz \
    > projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_126_R1.fq.gz

cat projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_112_R2.fq.gz \
    projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_113_R2.fq.gz \
    > projects/rbge/zedchen/barcoding/potamogeton_20250714/results/qc2_fastp/potamogeton_126_R2.fq.gz
}

function curate_easy353_assembly {

#cp -r /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/easy353_intermediate/phy3_easy353assembly/*\
#      /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/
#exclude pusillis6, pusillus2, alpinus3, alpinus5, crispus3, gramineus4,  praelongus1, perfoliatus6 and filiformis4/5, berchtoldi5 and berchtoldi6
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_106 -rf #pusillis6
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_46 -rf #pusillus2
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_63 -rf #alpinus3
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_78 -rf #alpinus5
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_79 -rf #crispus3
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_85 -rf #gramineus4
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_100 -rf #praelongus1
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_117 -rf #perfoliatus6
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_112 -rf #filiformis4
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_113 -rf #filiformis5
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_110 -rf #berchtoldi5
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_111 -rf #berchtoldi6
rm /mnt/shared/scratch/zchen/Barcoding_km/barcoding_potamogeton/results/phy3_easy353assembly/potamogeton_124 -rf #undetermined

}

#execution
function main {

#remove_combine
curate_easy353_assembly
}

main