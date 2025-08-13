#!/home/pholling/projects/rbge/pholling/env/snps/bin/python3
# 2025.08.12
# Zed ZChen@rbge.org.uk
# Input: vcf file, with GT:DP available (default)
# Output: csv file with genotype frequency for each species

from __future__ import division
import sys
import gzip
import copy
import argparse
import pandas as pd

def parse_vcf(vcf,output):
  #just in case the vcf file is gz
  if vcf.endswith(".gz"):
  	opener = gzip.open
  else:
  	opener = open
  #process each line
  snps,snps_info,n={},{},0
  for i in opener(vcf,'rt'):
    if i.startswith('#CHROM'):
      heading = i.strip().split()
      samples=heading[9:]
      samples=[s.split('/')[-1].replace('.bam','') for s in samples]
      print(f"processing {len(samples)} samples: {samples}")
    else:
      if not i.startswith('#'): #and n<=100; this is a good cap for testing just a small number of SNPs
        n+=1
        snps[f"SNP_{n}"]={s:{} for s in samples}
        #COLLECT SNP INFO AND OUTPUT IN A CSV
        snps_info[f"SNP_{n}"]={i:j for i,j in zip(heading[:9],i.strip().split()[:9])}
        #COLLECT SAMPLE INFO
        i = i.strip().split()[9:] #the sample starts from teh 9th column
        genotypes=[gt.split(':')[0] for gt in i]
        prob=[gt.split(':')[1] for gt in i]
        for s,gt in zip (samples,genotypes):
          snps[f"SNP_{n}"][s]['GT']=gt
        for s,pl in zip (samples,prob):
          snps[f"SNP_{n}"][s]['PL']=pl.split(',')
  #OUTPUT SNP INFO
  df=pd.DataFrame.from_dict(snps_info,orient='index') #turn to dataframe 
  df.to_csv(f'{output}/all_SNPs_summary.csv', index=True) 
  return(snps,snps_info)

#get samples vs species
def parse_names(names):
  id_sp={l.split(',')[0].split('/')[-1].replace('.bam',''):l.split(',')[1].strip('\n') for l in open(names)}
  return(id_sp)

#calcualte GT per species
def GT_per_species(id_sp,snps,snps_info,high,low,output):
  log=open(f"{output}/log.txt",'w')
  log.write(f"species specific SNP with a GT\n>= {high}% frequency in the present species,and \n<{low}% in absent species\n")
  sps=list(set(id_sp.values()))
  log.write(f"species: {sps}\n")
  sp_nb=len(sps)
  GTs=set([v2['GT'] for v in snps.values() for v2 in v.values()]) #extract all possible genotypes
  GTs={gt: f"GT_{i}" for i,gt in enumerate(GTs)}
  print(f"genotypes: {GTs}")
  log.write(f"genotypes: {GTs}\n")
  snp_GT={k:{} for k in snps.keys()}
  print(f"species included in this table: {sps}")
  for k,v in snps.items():
    for sp in sps:
      gts=[v[sample]['GT'] for sample,spp in id_sp.items() if sp == spp]
      snp_GT[k][sp]={gt:round(100*gts.count(gt)/len(gts),2) for gt in set(gts)} #now calculate the frequency
  GT_dict={gt:{} for gt in GTs}
  #look at each case of the GT: 0/0,0/1,1/1,etc.
  for gt,name in GTs.items():
    gt_freq={k:{species:freq.get(gt,0) for species,freq in v.items()} for k,v in snp_GT.items()} 
    log.write(f"process {len(gt_freq)} SNPs\n")
    GT_dict[gt]=gt_freq #the frequency of each genotype of each species
    df=pd.DataFrame.from_dict(gt_freq,orient='index') #turn to dataframe 
    df.to_csv(f'{output}/{name}_freq.csv', index=True) 
    print(df)
    #SELECT FOR SPECIES SPECIFIC GT OF THE SNP
    selected={}
    for k,v in gt_freq.items():
      freqs=[f for sp,f in v.items()]
      over=[1 for i in freqs if i >= high]
      belo=[1 for i in freqs if i <  low ]
      if sum(over)==1 and sum(over)+sum(belo)==sp_nb:
        selected[k]=v
        #add SNP info 
        selected[k]['genotype']=gt
        selected[k]={**selected[k],**snps_info[k]} #python3.5 or greater
    log.write(f"found {len(selected)} SNPs with species specific GT {gt}\n")
    df=pd.DataFrame.from_dict(selected,orient='index') #turn to dataframe 
    print(df)
    df.to_csv(f'{output}/{name}_freq_{low}_{high}.csv', index=True)  
  log.close()         
  return(GT_dict)
  
def main(vcf,names,high,low,output):
  #print(vcf,names,output)
  snps,snps_info=parse_vcf(vcf,output)
  id_sp=parse_names(names)
  GT_dict=GT_per_species(id_sp,snps,snps_info,int(high),int(low),output) #also a quick selection process
  
## files
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
	  description='calculate GT freq of each species in a vcf file', \
	  usage = 'calculate_snp.freq.py -v <vcf file> -n <ID_species.txt> -o <output_dir/prefix>')
  parser.add_argument('-v','--vcf', help='vcf file with GT and DP info',metavar='vcf')
  parser.add_argument('-n','--names', help='sample ID to species names corresponding file, csv: sample ID/file path, species',metavar='names')
  parser.add_argument('-i','--high', help='GT freq of species specific GT to be at least xx% in the present species',metavar='high')
  parser.add_argument('-l','--low', help='GT freq of species specific GT to be at most xx% in the absent species',metavar='low')
  parser.add_argument('-o','--output',help='path to output directory',metavar='output')
  options = parser.parse_args()
  
  main(options.vcf,options.names,options.high,options.low,options.output)
  
