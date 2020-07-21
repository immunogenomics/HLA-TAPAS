# MakeReference (v2.0)

## (1) Introduction

MakeReference builds a reference panel that can be subsequently used by [SNP2HLA](../SNP2HLA) for HLA imputation. It is a updated version of the original [MakeReference](http://software.broadinstitute.org/mpg/snp2hla/makereference_manual.html) module in SNP2HLA.
<br>

## (2) Input & output
You would need to files to get started with the panel making:
1. A file (in PLINK format) with all genotyped/sequenced genomic markers (\*.bed/bim/fam)
2. *HLA* alleles (after [NomenCleaner](../NomenCleaner)) for each individual included in 1 (\*.chped)

MakeReference will output a combined and **phased** HLA reference panel for imputing (in vcf format) and a maker file.

## (3) Usage example

MakeReference in HLA-TAPAS has to be implemented in the directory of main project folder. (i.e. 'HLA-TAPAS/' where 'HLA-TAPAS.py' script is.)

```
$ cd ../ 
# Change your current directory to the HLA-TAPAS main project folder.
```

```
$ python -m MakeReference \
    --variants MakeReference/example/HAPMAP_CEU \ # this is your marker file in the PLINK format
    --chped MakeReference/example/HAPMAP_CEU_HLA.imgt3320.4field.chped \ # this file includes all HLA types
    --hg 19 \
    --out MyRef/HAPMAP_CEU.REF.bglv4 \
    --dict-AA MakeReference/data/hg19/HLA_DICTIONARY_AA.hg19.imgt3320 \
    --dict-SNPS MakeReference/data/hg19/HLA_DICTIONARY_SNPS.hg19.imgt3320 \
    --phasing
```
