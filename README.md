
# VCF parse 

## Introduction

Generating a table annotating each variant in the VCF input file. 
Each variant is annotated with the following pieces of information:

1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API (API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that might be relevant.


## Requirements

This program requires the following modules/libraries (Python 2.7.12):

* numpy
* pandas
* requests
* json
* [vcf](https://pyvcf.readthedocs.io/en/latest/)
* [pyVEP](https://github.com/kantale/pyVEP)

## Documentation

Documentation is available for each of the functions inside the code (vcf-parse-challenge.ipynb). 

## Usage
```
python vcf-parse-challenge.py
```

## Author contact

* [Sara Movahedi](https://github.com/tabotaab) 
