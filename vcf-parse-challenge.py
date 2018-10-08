#!/usr/bin/env python
# coding: utf-8

'''
VCF parse 
==========

#  Introduction

Output a table annotating each variant in the VCF input file. Each variant must be annotated with the following pieces of information:
 
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API (API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that you feel might be relevant.
 
# Author contact

Sara Movahedi , https://github.com/tabotaab 
'''

# ## modules/libraries
# This module requires the following modules/libraries (Python 2.7.12):
# 
# * numpy
# * pandas
# * requests
# * json
# * [vcf](https://pyvcf.readthedocs.io/en/latest/)
# * [pyVEP](https://github.com/kantale/pyVEP)
# 


import vcf
import numpy as np
import pandas as pd
from pyVEP import VEP
import requests
import json


# ## Default variables



input_file = './Challenge_data.vcf'
output_file = './vcf-parse-output.csv'



info_type = {'del': "deletion",
              'ins': "insertion",
              'complex': "complex",
              'snp' : "Single Nucleotide Polymorphism",
              'mnp': "Multi Nucleotide Polymorphism",
              }
variation_type = {'indel' : "INDEL",
              'snp' : "SNP",
              }
variation_subtype = {'del': "deletion",
              'ins': "insertion",
              'ts': "transition",
              'tv': "transversion",
              'unknown': "unknown"}


# ## parse_vcf_record
# 
# A function to parse input vcf file records. 
# 
# **input:**
# * a vcf record = one line
# 
# **output:**
# * CHROM chromosome
# * POS position
# * REF reference base(s)
# * ALT comma separated list of alternate non-reference alleles called on at least one of the samples
# * var_type the type of variant (SNP or INDEL)
# * INFO TYPE comma separated list of variant types of alternate non-reference alleles
# * var_subtype variant subtype
# * AO comma separated list of alternate allele observation count
# * DP combined depth across samples
# * Percentage of reads supporting the variant versus those supporting reference reads



def parse_vcf_record(record):
    ref_reads, variant_reads, total_reads = record.INFO['RO'], record.INFO['AO'], record.INFO['DP']
                
    var_per = []
    for var_val in variant_reads:
        var_per.append(round((float(var_val)/float(total_reads))*100.0,2))
                    
    ref_per = round((float(ref_reads)/float(total_reads))*100.0,2)
                
    info_types = []
    for infotype in record.INFO['TYPE']:
        info_types.append(info_type[infotype])
                    
    tmp = (record.CHROM, record.POS, record.REF, str(record.ALT)[1:-1],
           variation_type[record.var_type],str(info_types)[1:-1],variation_subtype[record.var_subtype],
           str(variant_reads)[1:-1], total_reads, str(var_per)+"%|"+str(ref_per)+"%",)
                  
    return(tmp) 


# ## variant_effect
# 
# Here we use pyVEP to predict Variant Effect.
# 
# **input:**
# * a vcf record = one line
# 
# **output:**
# * comma separated list of 'most severe consequence' of alternate non-reference alleles
# 



def variant_effect(record): 
    
    var_effect = []
    
    for eachalt in record.ALT:
        mystr = str(record.CHROM)+" "+str(record.POS)+" . "+str(record.REF)+" "+str(eachalt)+" .  .  ."
        
        if len(record.REF)>1 and len(eachalt)>1:
            var_effect.append("complex")
        else:
            r = VEP(mystr, 'grch38')
            var_effect.append(str(r[0]['most_severe_consequence']))
           
    return((str(var_effect)[1:-1],))


# ## ExAC_info
# 
# Here we use Broad Institute ExAC Project API to gain more information over each allele.
# 
# **input:**
# * a vcf record = one line
# 
# **output:**
# * allele_count 
# * allele_num 
# * allele_freq 
# * num_homozygotes 
# * site_quality 
# * filter 
# * major_consequence
# 


def ExAC_info(record):
    chrom,pos,ref,alt = record.CHROM, record.POS, record.REF, record.ALT
    tmp = ("NA","NA","NA","NA","NA","NA","NA")
    
    if len(ref)==1 and len(alt)==1:
        exac_url = "http://exac.hms.harvard.edu/rest/variant/variant/"+chrom+'-'+str(pos)+'-'+ref+'-'+str(alt[0])
        exac_response = requests.get(exac_url)
        
        if (exac_response.status_code == 200) :
            exac_response_json = exac_response.json()
    
            allele_count = exac_response_json["allele_count"]
            allele_num = exac_response_json["allele_num"]
            allele_freq = exac_response_json["allele_freq"]
            num_homozygotes = exac_response_json["hom_count"]
            site_quality = exac_response_json["site_quality"]
            outFilter = exac_response_json["filter"]
    
            major_consequence={}
            vep_ann = exac_response_json["vep_annotations"]   
            for vep in vep_ann:
                major_consequence[vep["HGVSc"]]=vep["major_consequence"]
            
            tmp = (allele_count,allele_num,allele_freq,num_homozygotes,site_quality,str(outFilter),str(json.dumps(major_consequence))[1:-1])
        else:
            raise ValueError("ERROR: ExAC response status_code should be 200, but it is "+exac_response.status_code+" !")
        
    return(tmp)
    


# ## main function
# Here we read and parse the VCF input file line by line. Output results are written to a csv file. 


def main(vcf_file = input_file,out_file = output_file ):
    
    vcf_reader = vcf.Reader(open(vcf_file,'r'))

    f = open(out_file,"w+")
    f.write('chromosome'+'\t'+'position'+'\t'+'reference'+'\t'+'variant'+'\t'+'var_type'+'\t'+'var_infotype'+'\t'+
            'var_subtype'+'\t'+'var_count'+'\t'+'read_depth'+'\t'+'var%|ref%'+'\t'+'var_effect'+'\t'+
            'ExAC_allele_count'+'\t'+'ExAC_allele_num'+'\t'+'ExAC_allele_freq'+'\t'+'ExAC_homozygotes_num'+'\t'+
            'ExAC_site_quality'+'\t'+'ExAC_filter'+'\t'+'ExAC_HGVSc:major_consequence'+'\n')
    
    for record in vcf_reader:
        try:
                vcf_out = parse_vcf_record(record) # VCF information
                eff_out = variant_effect(record)   # Variant effect results
                
                try:
                    ExAC_out = ExAC_info(record)   # ExAC API output
                except KeyError:
                    ExAC_out = ("NA","NA","NA","NA","NA","NA","NA")
                
                out_line = vcf_out+eff_out+ExAC_out
                f.write('\t'.join(str(i) for i in out_line)+'\n')
        except KeyError:
            print('WARNING: missing count field(s) in record %s:%d' % (record.CHROM, record.POS))

    f.close()
    
    return()



if __name__== "__main__":
  main()

