#!/home/zchen/apps/env/easy353_2/bin/python3

import subprocess as sbp
import sys
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#list of functions
func_list=[
  'open_fasta',
  'collect_all_ref',
  'extract_same_chrom',
  'check_homology'
]

#check each taxa
def open_fasta(fafile):
  fadict=SeqIO.to_dict(SeqIO.parse(fafile,'fasta'))
  return(fadict)

def collect_all_ref(indir):
  falist=sbp.check_output(f"ls {indir}/*fasta",shell=True).decode().strip('\n').split('\n')
  fadicts={fa.split('/')[-1].replace('fasta',''):open_fasta(fa) for fa in falist}
  return fadicts

#if labelled the same chromosome, put them into the same fasta file for subsequent blasting confirmation
def extract_same_chrom(indir,outdir):
  fadicts=collect_all_ref(indir)
  for i in range(1,20):
    print(f'chromosome: {i}')
    records=[]
    for sp,fadict in fadicts.items():
      for record in fadict.values():
        if record.description.endswith(f'chromosome: {i}'):
          print(sp,record.description)
          records.append(record)
    SeqIO.write(records,f'{outdir}/chromosome_{i}.fasta','fasta')

def check_homology(indir):
  for i in range(1,20):
    print(f"blasting chromosome {i} of the reference species")
    #sbp.call(f"./blast_two_files.sh {indir}/chromosome_{i}.fasta {indir}/chromosome_{i}.fasta 6", shell=True)
    sbp.call(f"makeblastdb -in {indir}/chromosome_{i}.fasta -dbtype nucl -out tmp_chromosome_{i}", shell=True)
    sbp.call(f"blastn -db tmp_chromosome_{i} -query {indir}/chromosome_{i}.fasta -outfmt 6 > {indir}/blast_chromosome_{i}_out.txt", shell=True)
    sbp.call(f"rm tmp_chromosome_{i}* -f",shell=True)
    print(f"Completed blasting chromosome {i} of the reference species")

def split_species_chrom(fafile,outdir):
  fadict=SeqIO.to_dict(SeqIO.parse(fafile,'fasta')) 
  for i in range(1,20):
    records=[]
    for record in fadict.values():
      if record.description.endswith(f'chromosome: {i}'):
        print(record.description)
        records.append(record)
    SeqIO.write(records,f'{outdir}/Sh_chromosome_{i}.fasta','fasta')

def check_length(indir):
  fadicts=collect_all_ref(indir)
  for sp,fadict in fadicts.items():
    for record in fadict.values():
      print(sp,len(record.seq),record.description,)
#EXECUTE FUNCTIONS
def main(indir,outdir):
  fafile='/mnt/shared/scratch/zchen/Barcoding_km/SNP_potamogeton/refs/S_herbacea.chromosome.fasta' #use herbacea as the bait
  #extract_same_chrom(indir,outdir)
  #check_homology(indir)
  #split_species_chrom(fafile,outdir)
  check_length(indir)
##################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='check references', \
	  usage = 'process_ref.py <indir> <outdir>')
  parser.add_argument('indir', help='directory of the references', metavar='indir')
  parser.add_argument('outdir', help='directory to output the chromosome files', metavar='indir')
  options = parser.parse_args()
  
  main(options.indir,options.outdir)