#!/usr/bin/env python3

from os import system
import subprocess as sbp
import argparse
from Bio import SeqIO
from Bio import Align
from Bio.SeqRecord import SeqRecord

###############################################
#standard input: ./detect_conserved_SSR.py 
###############################################

def read_SSR_to_dict(misa,prefix):
  with open(misa) as handle:
    keys=[x for x in next(handle).strip('\n').split('\t') if x !='']
    misa_dict={f'{prefix}_SSR_{i}': {i:j for i,j in zip(keys,[x for x in line.strip('\n').split('\t') if x !=''])}for i,line in enumerate(handle)}
    misa_dict={k:{i:int(j) if i in ['size', 'start', 'end'] else j for i,j in v.items() } for k,v in misa_dict.items()} #make sure int entries are actually numbers
  #print(misa_dict)
  return(misa_dict)  

def run_blastn(contigs1,misa1,prefix1,contigs2,misa2,prefix2,flank):
  #use subprocess to fun a few scripts in the commandline
  sbp.call(' '.join(['./extract_flanked.py',contigs1,misa1,prefix1,flank]),shell=True) #extract SSRs from the first group
  sbp.call(' '.join(['./extract_flanked.py',contigs2,misa2,prefix2,flank]),shell=True) #extract SSRs from the second group
  #(added commandline_blastn.sh path to .bashrc file, can be reused later)
  sbp.call(' '.join(['./command_line_blast.sh',f"{prefix1}_SSR_flk_{flank}.fasta",f"{prefix2}_SSR_flk_{flank}.fasta",'6']),shell=True) #blast the two files 
  #prefix1: db, prefix2: query
  with open('blast_out.txt') as handle:
    keys=['qseqid', 'sseqid' ,'pident', 'length' ,'mismatch' ,'gapopen', 'qstart' ,'qend', 'sstart', 'send', 'evalue', 'bitscore']
    blast_dict={f'pair_{i}': {i:j for i,j in zip(keys,[x for x in line.strip('\n').split('\t') if x !=''])}for i,line in enumerate(handle)}
    blast_dict={k:{i:float(j) if i in ['pident', 'length' ,'mismatch' ,'gapopen', 'qstart' ,'qend', 'sstart', 'send', 'evalue', 'bitscore'] 
                              else j for i,j in v.items() } 
                              for k,v in blast_dict.items()} #make sure int entries are actually numbers
  return(blast_dict) 

def add_seqs(misa_dict,prefix,flank):
  seqs=SeqIO.to_dict(SeqIO.parse(f"{prefix}_SSR_flk_{flank}.fasta",'fasta'))
  misa_dict={k:v for k,v in misa_dict.items() if k in seqs}
  for k in seqs:
    misa_dict[k]['seq']=seqs[k].seq
  return misa_dict

def remove_number_in_string(string): 
#compare only the repeat themselves
  string=[i for i in string if i not in ['0','1','2','3','4','5','6','7','8','9']]
  return ''.join(string)

def group_SSR_pairs(blast_dict,misa_dict1,misa_dict2):
#look for different SSRs in a conserved flanking region
  identical_SSRs,diff_SSRs=[],[] #identical+flank, same unit but different number of repeats, different units
  for k,v in blast_dict.items():
    if v['pident']==100:
      identical_SSRs.append((k,blast_dict[k]['sseqid'],blast_dict[k]['qseqid'])) #in order of misa_dict1,misa_dict2
      #print(f"Discard {k}: 100% PID") 
    #prefix1: db, prefix2: query
    else:
      ssr1=str(misa_dict1[blast_dict[k]['sseqid']]['seq'][150:-150])
      ssr2=str(misa_dict2[blast_dict[k]['qseqid']]['seq'][150:-150])
      if ssr1==ssr2:
        #print(ssr1,ssr2)
        identical_SSRs.append((k,blast_dict[k]['sseqid'],blast_dict[k]['qseqid']))
      else:
        diff_SSRs.append((k,blast_dict[k]['sseqid'],blast_dict[k]['qseqid']))
      #print(f"{misa_dict2[blast_dict[k]['qseqid']]['SSR']}!={misa_dict1[blast_dict[k]['sseqid']]['SSR']}")
  stats=['\n###############################\n','Summary of SSR blastn pairs',
          f'Total: {len(blast_dict)}',
          f'Identical Sequence: {len(identical_SSRs)}',
          f'Same SSR unit different number of repeats: {len(diff_SSRs)}',
          '\n###############################\n']
  print('\n'.join(stats))
  return {'identical':identical_SSRs, 'diff_SSRs':diff_SSRs}

def filter_SSR_pairs(blast_dict,misa_dict1,misa_dict2,flank): #group is a list in the dict 'groups'
  count1,count2,count3,cutoff=0,0,0,flank*2-50
  tmp={k:v for k,v in blast_dict.items() if v['length']>=cutoff} #must align for over certain bp
  f1=open(f'compare_SSR_identical_{flank}.txt','w')
  f2=open(f'compare_SSR_same_unit_{flank}.txt','w') #to collect output
  f3=open(f'compare_SSR_diff_unit_{flank}.txt','w')
  for k,v in tmp.items():
    ssr1,ssr2=v['sseqid'],v['qseqid']
    #reorient the sequences
    if v['sstart']> v['send'] and v['qstart'] > v['qend']: #reorient
      seq1=misa_dict1[blast_dict[k]['sseqid']]['seq'].reverse_complement() #if mapped in minus, reverse complement
      seq2=misa_dict2[blast_dict[k]['qseqid']]['seq'].reverse_complement() #if mapped in minus, reverse complement
    elif v['sstart']>v['send'] and v['qstart']<v['qend']: 
      seq1=misa_dict1[blast_dict[k]['sseqid']]['seq'].reverse_complement() #if mapped in minus, reverse complement
      seq2=misa_dict2[blast_dict[k]['qseqid']]['seq']
    elif v['sstart']<v['send'] and v['qstart']>v['qend']:
      seq1=misa_dict1[blast_dict[k]['sseqid']]['seq']
      seq2=misa_dict2[blast_dict[k]['qseqid']]['seq'].reverse_complement() #if mapped in minus, reverse complement
    elif v['sstart']<v['send'] and v['qstart']<v['qend']:
      seq1=misa_dict1[blast_dict[k]['sseqid']]['seq']
      seq2=misa_dict2[blast_dict[k]['qseqid']]['seq']
  #finished reorienting: flip one, flip both, not flip
    seqs1='\t'.join([ssr1,str(seq1[:flank]),str(seq1[flank:-flank]),str(seq1[-flank:])])
    seqs2='\t'.join([ssr2,str(seq2[:flank]),str(seq2[flank:-flank]),str(seq2[-flank:])])
    if str(seq1[flank:-flank])==str(seq2[flank:-flank]):
      f1.write('\n'.join([k,seqs1,seqs2,'\n']))
      count1+=1
    else:
      
      if set(list(str(seq1[flank:-flank])))==set(list(str(seq2[flank:-flank]))): #if the repeats only differ in the number of repeats
        f2.write('\n'.join([k,seqs1,seqs2,'\n']))
        count2+=1
      else:
        f3.write('\n'.join([k,seqs1,seqs2,'\n']))
        count3+=1
  print('\n'.join(['\n###############################\n',
        f"Total SSR pairs with >={cutoff}bp alignment: {len(tmp)}",
        f"Number of identical SSR pairs after reorienting: {count1}",
        f"Number of SSR pairs with identical units but different number of repeats: {count2}",
        f"Number of SSR pairs with different units : {count3}",
        '\n###############################\n']))
  f1.close()
  f2.close()
  f3.close()

def main(contigs1,misa1,prefix1,contigs2,misa2,prefix2,flank=150):
  misa_dict1=read_SSR_to_dict(misa1,prefix1)
  misa_dict2=read_SSR_to_dict(misa2,prefix2)
  blast_dict=run_blastn(contigs1,misa1,prefix1,contigs2,misa2,prefix2,flank) #get the blastn results
  #add sequences to selected SSRs
  misa_dict1=add_seqs(misa_dict1,prefix1,flank)
  misa_dict2=add_seqs(misa_dict2,prefix2,flank)
  #groups=group_SSR_pairs(blast_dict,misa_dict1,misa_dict2) #group the blast pairs
  #process some groups
  filter_SSR_pairs(blast_dict,misa_dict1,misa_dict2,int(flank)) #aligned for over 250 bp
  #print(blast_dict)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='detect conserved SSRs')
  parser.add_argument('contigs1', help='Name of the genome assembly, filtered.fasta', metavar='contigs1')
  parser.add_argument('misa1', help='Name of the misa output', metavar='misa1')
  parser.add_argument('prefix1', help='Prefix of the SSR output', metavar='prefix1')
  parser.add_argument('contigs2', help='Name of the genome assembly, filtered.fasta', metavar='contigs2')
  parser.add_argument('misa2', help='Name of the misa output', metavar='misa2')
  parser.add_argument('prefix2', help='Prefix of the SSR output', metavar='prefix2')
  parser.add_argument('flank', help='Required length of flanking region of the SSR', metavar='flank')
  options = parser.parse_args()
  
  main(options.contigs1, options.misa1,options.prefix1, options.contigs2, options.misa2,options.prefix2,options.flank)