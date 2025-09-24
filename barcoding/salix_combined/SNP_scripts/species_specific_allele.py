#!/home/zchen/apps/env/easy353/bin/python3
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
import subprocess as sbp

#get samples vs species
def parse_names(names):
  id_sp={l.split(',')[0].split('/')[-1].replace('.bam',''):l.split(',')[1].strip('\n') for l in open(names)}
  return(id_sp)

def parse_single_record(l,samples,heading):
  snps={s:{} for s in samples}
  #COLLECT SNP INFO AND OUTPUT IN A CSV
  snps_info={i:j for i,j in zip(heading[:9],l.strip().split()[:9])}
  #COLLECT SAMPLE INFO
  l = l.strip().split()[9:] #the sample starts from teh 9th column
  genotypes=[gt.split(':')[0] for gt in l]
  for s,gt in zip (samples,genotypes):
    snps[s]['GT']=gt
  return(snps,snps_info)

#calcualte GT per species
def species_specific_alleles(k,id_sp,sps,samples,snps,snps_info,hqallele,hqalGT,high):
  snp_al,snp_GT={},{}
  for sp in sps: #process each species
    gts=[snps[sample]['GT'] for sample,spp in id_sp.items() if sp == spp and sample in samples]
    alleles=[]
    for gt in gts:
      alleles+=[int(i) for i in gt.split('/') if i != '.'] #remove the uncertain ones and collect all alleles
    #calculate allele frequency of the sample
    snp_al[sp]={al:round(alleles.count(al)/len(alleles)*100) for al in set(alleles)} #store the value
    #calculate GT frequency
    gts=[snps[sample]['GT'] for sample,spp in id_sp.items() if sp == spp and sample in samples]
    snp_GT[sp]={gt:round(100*gts.count(gt)/len(gts),2) for gt in set(gts)} #the dictionary with all GT frequencies for all SNPs
  tmp1={f"{species}_{al}": v.get(al,0) for species,v in snp_al.items() for al in [0,1]}
  tmp2={f"{species}_{gt}": v.get(gt,0) for species,v in snp_GT.items() for gt in ['0/0','0/1','1/1']}
  #print(snp_GT[k])
  #filter for SNPs with HQ species specific allele
  for al,typ in {0:'REF',1:'ALT'}.items(): #assume two alleles
    freqs={sp:snp_al[sp].get(al,0) for sp in snp_al}
    present=[1 for f in freqs.values() if f > 0 ]
    if sum(present)==1 : #the allele is present in only one species
      ssallele,target={},[s for s,f in freqs.items() if f >0][0]
      ssallele['target sp']=target
      ssallele['target allele']=snps_info[typ]
      ssallele['target freq%']=freqs[target]
      ssallele={**ssallele,**snps_info} #add snps info
      ssallele={**ssallele,**tmp1} #add allele freq for each species
      #print(k,al, sum(freqs),high) #make the first dictionary for collecting values
      if sum(freqs.values())>=high:
        print(k,ssallele['target sp'],ssallele['target freq%'])
        hqallele[k]=ssallele #record if hq
        hqalGT[k]={**ssallele,**tmp2} #add GT freq info
        #hqal.append(','.join(ssallele.values()))
        #hqal.append('\n')
        continue #no need to check the other allele
      else:
        continue
      

#DONE WITH THE TRAINING DATASET
#ONTO THE QUERY SET
#EXTRACT GOOD SNPS
def extract_SNPs(snps,snps_info,snp_out):
  #reset the keys for extraction
  print('number of SNPs before renaming keys',len(snps),len(snps_info))
  snps={f"{v['#CHROM']}_{v['POS']}":snps[k] for k,v in snps_info.items()}
  snps_info={f"{v['#CHROM']}_{v['POS']}":v for k,v in snps_info.items()}
  print('number of SNPs after renaming keys',len(snps),len(snps_info))
  #extract
  print('number of species specific SNPs',len(snp_out)) 
  snps={k:snps.get(f"{v['#CHROM']}_{v['POS']}",'NA') for k,v in snp_out.items()}
  snps_info={k:snps_info.get(f"{v['#CHROM']}_{v['POS']}",'NA') for k,v in snp_out.items()}
  print('number of SNPs after extraction',len(snps),len(snps_info))
  return(snps,snps_info)

#SHOW SAMPLE GT
def specific_SNPs_GT(samples,snps,snps_info,snp_out,output):
  GTs={}
  for k,v in snp_out.items():
    if snps_info[k]=='NA': #if the snp is not called in those samples
      print(f'{k} not called in all samples: identical to reference')
    else:
      #print(snps_info[k],v)
      if snps_info[k]['ALT']==v['ALT']:
        GTs[k]={s:snps[k][s]['GT'] for s in samples}
      else:
        for s in snps[k]:
          if snps[k][s]['GT']=='0/0':
            GTs[k][s]='0/0'
          else:
            GTs[k]='novo'
  df=pd.DataFrame.from_dict(GTs,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_SNP_sample_GT.csv', index=True)  
  #convert GTs to list of alleles
  allele_table={0:'REF',1:'ALT'}
  alleles={snp:{} for snp in GTs}
  for snp,v in GTs.items():
    alleles[snp]={sample:[snps_info[snp][allele_table[int(i)]] for i in gt.split('/')] for sample,gt in v.items()}
    alleles[snp]['target sp']=snp_out[snp]['target sp']
    alleles[snp]['target allele']=snp_out[snp]['target allele']
  df=pd.DataFrame.from_dict(alleles,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_SNP_sample_allele.csv', index=True)  
  return(GTs,alleles)

#TRANSLATE SAMPLE GT TO SPECIES
def GT_to_species(alleles,snp_out,output):
  species={}
  for k,v in alleles.items():
    species[k]={sample:snp_out[k]['target sp'] if snp_out[k]['target allele'] in al else 'NA' for sample,al in v.items()}
    species[k]['target sp']=snp_out[k]['target sp']
    species[k]['target allele']=snp_out[k]['target allele']
  df=pd.DataFrame.from_dict(species,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/specific_allele_sample_species.csv', index=True)  
  return(species)

#
#
#
#PARENT FUNCTION 1 : mem efficient loop version
#using samples with known species, identify species specific SNPs
def specifi_SNPs(vcf,names,high,low,output):
  id_sp=parse_names(names) #sample names
  ssallele,hqallele,hqalGT={},{},{} #data collector
  #hqal=open(f'{output}/hq_specific_allele_freq.csv','w')
  #hqalgt=open(f'{output}/hq_specific_allele_GT_freq.csv','w')
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
      sps=list(set([sp for sample,sp in id_sp.items() if sample in samples])) #only consider species present in the sample. the sample-species list can be longer/contain more species than the one analysed
      sp_nb=len(sps)
    else:
      if not i.startswith('#'): 
        n+=1
        snps,snps_info=parse_single_record(i,samples,heading) #parse one line
        #calculate allele freq and collect hq SNPs
        species_specific_alleles(f"SNP_{n}",id_sp,sps,samples,snps,snps_info,hqallele,hqalGT,int(high))
  #output all HQ records
  #hqal.close()
  #hqalgt.close()
  df=pd.DataFrame.from_dict(hqallele,orient='index') #turn to dataframe 
  print(df)
  df.to_csv(f'{output}/hq_specific_allele_freq.csv', index=True)  
  df=pd.DataFrame.from_dict(hqalGT,orient='index') #output allele as well as GT freq
  print(df)
  df.to_csv(f'{output}/hq_specific_allele_GT_freq.csv', index=True)  



#PARENT FUNCTION 2 #haven't converted to loop version
#for samples without known species, find species/potential parents of hybrids based on species specific SNPs
def find_species(vcf,names,high,low,output):
  #read species specific SNPs from the previous step
  df=pd.read_csv(names,index_col=0)
  print(df)
  snp_out=df.to_dict(orient='index')
  #COLLECT SNPS
  snps,snps_info,samples=parse_vcf(vcf,output)
  #EXTRACT SPECIES SPECIFIC SNPS (IF EXIST)
  snps,snps_info=extract_SNPs(snps,snps_info,snp_out)
  GTs,alleles=specific_SNPs_GT(samples,snps,snps_info,snp_out,output)
  species=GT_to_species(alleles,snp_out,output)

#now making the execution function  
def main():
  count=0
  #parser = argparse.ArgumentParser(description='rename files of a certain type within the directory to use slurm array later',  usage = 'rename_files.py -i <indir> -o <outdir> --prefix <pre>')
  parser = argparse.ArgumentParser(description="Python script with multiple executable functions")
  subparsers = parser.add_subparsers(dest='command', help='Available commands')
  
  # Common arguments that will be reused
  common_args = argparse.ArgumentParser(add_help=False)
  common_args.add_argument('-v','--vcf', help='vcf file with GT and DP info',metavar='vcf')
  common_args.add_argument('-n','--names', help='sample ID to species names corresponding file, csv: sample ID/file path, species',metavar='names')
  common_args.add_argument('-i','--high', help='GT freq of species specific GT to be at least xx% in the present species',metavar='high')
  common_args.add_argument('-l','--low', help='GT freq of species specific GT to be at most xx% in the absent species',metavar='low')
  common_args.add_argument('-o','--output',help='path to output directory',metavar='output')
  
  #find species specific SNPs
  func_parser=subparsers.add_parser('specifi_SNPs', parents=[common_args],help='find and output species specific SNPs', 
                                    usage = './calculate_snp.freq.py main_parent -v <vcf file> -n <ID_species.txt> -i <high bound> -l <low bound> -o <output_dir/prefix> ')
  func_parser.set_defaults(func=specifi_SNPs)
  
  #assign species to unclassified samples based on SNPs
  func_parser=subparsers.add_parser('find_species', parents=[common_args],help='find and output hits on identified species specific SNPs', 
                                    usage = './calculate_snp.freq.py find_species -v <vcf file> -n <species specific SNPs.csv> -i <high bound> -l <low bound> -o <output_dir/prefix> ')
  func_parser.set_defaults(func=find_species)
  
  #parse arguments
  args = parser.parse_args()

  if not args.command:
    parser.print_help()
    sys.exit(1)
  
  # Prepare common kwargs for the function call
  kwargs = {'vcf': args.vcf,
            'names': args.names,
            'high':args.high,
            'low':args.low,
            'output':args.output}
  
  # Add function-specific arguments
  #if args.command == 'rename_fqgz':
  #  kwargs['prefix'] = args.prefix
  #elif args.command == 'rename_fasta':
  #  kwargs['outfile'] = args.outfile
  #elif args.command == 'rename_contig':
  #  kwargs['infile'] = args.infile
  #  kwargs['fcsv'] = args.fcsv
  
  # Call the appropriate function
  args.func(**kwargs)
    

if __name__ == '__main__': 
  main()