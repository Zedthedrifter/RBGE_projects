#!/home/zchen/projects/rbge/zedchen/env/easy353/bin/python3

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
  'process_tree',
  'process_treefiles',
  'is_monophyletic',
  'gene_summary'
]

#check each taxa
def is_monophyletic(tree, prefix):
  target_tips = [tip.name for tip in tree.get_terminals() if tip.name.startswith(prefix)]
  if not target_tips:
    return False
    # Find the MRCA of these tips
  try:
    mrca = tree.common_ancestor(target_tips)
  except:
    # This can happen if some tips are not found (shouldn't happen with our filtering)
    return False
  clade_tips = [tip.name for tip in mrca.get_terminals()]
  return all(tip.startswith(prefix) for tip in clade_tips) #check if all members start with the same prefix (if yes, monophyletic)

#check each tree
def process_tree(treefile,mono_dict):
  try:
    tree = Phylo.read(treefile, 'newick')
  except:
    print(f"Error reading tree file: {treefile}")
    return 0
  gene=treefile.split('/')[-1].replace('.fasta.treefile','')
  tip_names = [tip.name for tip in tree.get_terminals() if tip.name]
  pres=set([''.join([i for i in name if i not in ['0','1','2','3','4','5','6','7','8','9']]) for name in tip_names])
  for prefix in pres:
    if is_monophyletic(tree, prefix):
      mono_dict[gene][prefix]=1
    else:
      mono_dict[gene][prefix]=0

#iterate the whole directory
def process_treefiles(dirc,outdir):
  treefiles=sbp.check_output(f"ls {dirc}/*/*treefile",shell=True).decode().strip('\n').split('\n')
  mono_dict={treefile.split('/')[-1].replace('.fasta.treefile',''):{} for treefile in treefiles}
  for treefile in treefiles:
    process_tree(treefile,mono_dict) #add value to key
  df=pd.DataFrame.from_dict(mono_dict, orient='index')
  df=df.fillna(0)
  df.to_csv(f'{outdir}/gene_resolution_local.csv', index=True) 
  print(df)

#summary the power to resolve monophyletic group from each gene
def gene_summary(dirc,outdir):
  df=pd.read_csv(f'{outdir}/gene_resolution_local.csv',index_col=0)
  names=pd.read_csv(f'{outdir}/sample_names.csv',index_col=0)
  #row sum: the number of taxa each gene can resolve
  df['total_taxa']=df.sum(axis=1, skipna=True)
  df = df.sort_values('total_taxa', ascending=False)
  print(df['total_taxa'].value_counts()) #a summary of genes with different resolving capacity
  #select for genes with high resolving capacity
  for cut in [1]: #can only resolve one taxa
    try:
      sbp.call(f"mkdir {outdir}/gene_{cut}taxa",shell=True)
    except:
      print(f"{outdir}/gene_{cut}taxa already exist")
    #concatenate the genes
    concat={k:'' for k in names.index } #reset concatenation collector
    count=0
    #CONCATENATION
    for gene in df[df['total_taxa']>=cut].index:#select for genes with at least xx resolving power
      print(gene)
      contigs=SeqIO.to_dict(SeqIO.parse(f"{dirc}/{gene}/{gene}.fasta",'fasta')) #get from the phylogeny directory
      genlen=[len(v.seq)for k,v in contigs.items()][0] #len of the aligned gene
      for k,v in concat.items(): #add to each species
        try:
          concat[k]=v+contigs[k].seq
        except:
          concat[k]=v+Seq('-'*genlen)#if a certain sample doesn't have that gene, just fill with gap of the same length
      count+=1
    print(f"added {count} orthogroups")
    #print(concat)
    records=[SeqRecord(v,id=k,name=k,description=f"genes that can resolve {cut} taxa") for k,v in concat.items() if len(set(str(v))) != 1] #remove all gap cases
    lens=[len(v.seq) for v in records]
    print('concatenated contig lenghts:\n',lens)
    SeqIO.write(records,f"{outdir}/gene_{cut}taxa/gene_{cut}taxa.fasta",'fasta')
    sbp.call(f'iqtree -s {outdir}/gene_{cut}taxa/gene_{cut}taxa.fasta -bb 1000 -redo -safe',shell=True)

#select for a set of genes that resolve the most monophyletic taxa 
def gene_select(dirc,outdir,cvg):
  df=pd.read_csv(f'{outdir}/gene_resolution_local.csv',index_col=0)
  names=pd.read_csv(f'{outdir}/sample_names.csv',index_col=0)

#FUNCTIONS
def main(dirc,outdir):
  #process_treefiles(dirc,outdir) #only need to run once
  #gene_summary(dirc,outdir)  
  #row sum: the number of taxa each gene can resolve
  gene_select(dirc,outdir,cvg)

##################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='output a dataframe for the number of monophyletic groups that can be resolved', \
	  usage = 'count_mono.py <dir> <outdir>')
  parser.add_argument('dirc', help='parent directory to 353 gene phylogenies', metavar='dirc')
  parser.add_argument('outdir', help='output directory', metavar='outdir')
  options = parser.parse_args()
  
  main(options.dirc,options.outdir)