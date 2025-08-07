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
  'open_tree' #rename fqgz files to use array
]

#FUNCTIONS
def open_tree(treefile):
  tree=open(treefile).read()
  print(tree)
  return(tree)

#filter out bad taxa and bad seq
def swap_names(tree,t1,t2):
  print('swapping')
  tree=tree.replace(t1,'hold').replace(t2,t1).replace('hold',t2)
  return(tree)
##################################################
def main(treefile):
  tree=open_tree(treefile)
  tree=swap_names(tree,'P_gramineus_4','P_crispus_3')  #PG1 vs PC1
  tree=swap_names(tree,'S_filiformis_4','P_perfoliatus_6')  #SEPA_14a vs SEPA_4
  tree=swap_names(tree,'S_filiformis_5','P_perfoliatus_6b')  #SEPA_14b vs SEPA_4
  tree=swap_names(tree,'P_alpinus_5','P_praelongus_1')  #PAL1 & PPR1
  tree=swap_names(tree,'P_alpinus_3','P_lucens_6')
  print(tree)
  f=open(treefile.replace('treefile','corrected.treefile'),'w')
  f.write(tree)
  f.close()

##################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='remove the taxa that are not good for phylogeny construction', \
	  usage = 'tree_modifier.py <treefile>')
  parser.add_argument('treefile', help='treefile to be modified', metavar='treefile')
  options = parser.parse_args()
  
  main(options.treefile)