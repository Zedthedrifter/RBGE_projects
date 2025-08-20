#!/home/zchen/apps/env/easy353_2/bin/python3

import subprocess as sbp
import sys
import argparse
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#list of functions
func_list=[
  'remove_taxa' #rename fqgz files to use array
]

#FUNCTIONS

#filter out bad taxa and bad seq
def remove_taxa(fasta,bad_smp):
  for seq in SeqIO.parse(fasta,'fasta'): #parse can only be used as iterator once. you forgot that
    print(seq.id)
    if seq.id in bad_smp:
      print(seq.id)
  #the actual output
  seqs=SeqIO.parse(fasta,'fasta')
  records=[SeqRecord(seq.seq,id=seq.id.replace(' ','').replace('.','_'),name='',description='') 
  for seq in seqs 
  if seq.id.replace(' ','').replace('.','_') not in bad_smp and len(set(str(seq.seq))) != 1]
  SeqIO.write(records,fasta,'fasta')
  

##################################################
def main(fasta):  
  bad_smp=['P_pusillus_6','P_pusillus_2','P_alpinus_3','P_alpinus_5','P_crispus_3',
           'P_gramineus_4','P_praelongus_1','P_perfoliatus_6','S_filiformis_4',
           'S_filiformis_5','P_berchtoldii_5','P_berchtoldii_6','undetermined']
  remove_taxa(fasta,bad_smp)
    

##################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='remove the taxa that are not good for phylogeny construction', \
	  usage = 'remove_taxa.py <fasta> <taxa>')
  parser.add_argument('fasta', help='fasta file to be filtered', metavar='fasta')
  options = parser.parse_args()
  
  main(options.fasta)
  
#just some random functions not used in this script
def rename_ref(fasta):
  records=[]
  for c in SeqIO.parse(fasta,'fasta'):
    heading=c.id.split('-')
    heading=f"seq{heading[1]}-{heading[0]}"
    r=SeqRecord(c.seq,id=heading,description='',name='')
    records.append(r)
  SeqIO.write(records,'out.fasta','fasta')