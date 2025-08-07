#!/home/zchen/projects/rbge/zedchen/env/easy353/bin/python3

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
def move_reindex(indir,outdir,add):
  r1=sbp.check_output(f"ls {indir}/*1.fq.gz",shell=True).decode().strip('\n').split('\n')
  r2=sbp.check_output(f"ls {indir}/*2.fq.gz",shell=True).decode().strip('\n').split('\n')
  #step1: get new names by reindexing
  index1=[r.split('/')[-1].replace('salix_','').replace('_1.fq.gz','').replace('R','') for r in r1]
  index2=[r.split('/')[-1].replace('salix_','').replace('_2.fq.gz','').replace('R','') for r in r2]
  if set(index1)==set(index2):
    print('indexing is coordinated')
    new1={i:i.split('/')[-1].replace(j,str(int(j)+add)).replace('R','') for i,j in zip(r1,index1)}
    new2={i:i.split('/')[-1].replace(j,str(int(j)+add)).replace('R','') for i,j in zip(r2,index2)}
    #print(new1,'/n',new2)
  #step2: move and rename the files
    for k,v in new1.items():
      sbp.call(f"mv {k} {outdir}/{v}",shell=True)
    for k,v in new2.items():
      sbp.call(f"mv {k} {outdir}/{v}",shell=True)
  else:
    print(index1,index2)

def renaming_1_to_R1(indir):
  r1=sbp.check_output(f"ls {indir}/*1.fq.gz",shell=True).decode().strip('\n').split('\n')
  r2=sbp.check_output(f"ls {indir}/*2.fq.gz",shell=True).decode().strip('\n').split('\n')
  new1={i:i.replace('_1.fq.gz','_R1.fq.gz') for i in r1}
  new2={i:i.replace('_2.fq.gz','_R2.fq.gz') for i in r2}
  #print(new1)
  #print(new2)
  for k,v in new1.items():
    sbp.call(f"mv {k} {v}",shell=True)
  for k,v in new2.items():
    sbp.call(f"mv {k} {v}",shell=True)

def job_test():
  print('running job')
##################################################
def main():
  indir1='$SCRATCH/Barcoding_km/barcoding_salix_28/results/qc2_fastp'
  indir2='$SCRATCH/Barcoding_km/barcoding_salix_13/results/qc2_fastp'
  outdir='$SCRATCH/Barcoding_km/barcoding_combined_salix/results/00_reads'
  #move_reindex(indir1,outdir,47)
  #move_reindex(indir2,outdir,75)
  renaming_1_to_R1(outdir)
  job_test()
main()