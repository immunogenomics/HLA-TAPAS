# NomenCleaner

## (1) Introduction

`NomenCleaner` is a module to help researchers to solve challenges related to HLA nomenclature. Two major challenges would be (1) Conversion between nomenclatures, and (2) Field checking.



1. Conversion

NomenCleaner converts a given set of HLA alleles which is in arbitrary nomenclature to the ones in the updated nomenclature. **The conversion from the updated to old nomenclature is NOT available**.

<br>

2. Field checking

**Often, the field separator is removed though given HLA allele conforms to the updated nomenclature**. This barely causes a problem in most cases, however, there is definitely some case that may confuse researchers. For example, if the HLA allele **DPB1\*101101** is given, then some researchers may be perplexed as to which one is right among <U>DBP1\*10:11:01</U>, <U>DPB1\*101:101</U> and <U>DPB1\*1011:01</U>. NomenCleaner searches the `HAT(HLA Allele Table)` file, which contains whole HLA allele name information of the IMGT database of a specific version, and determines that **DPB1*1011:01 is the only valid solution**. (cf. DBP1*1011:01 is in the IMGT database of the version 3.37.0.)

<br>

## (2) Input & output
NomenCleaner basically takes `PED` file(e.g. '\*.ped') and `HAT(HLA Allele Table)` file (e.g. '\*.hat') as input, then performs its jobs as required by the user. Finally the NomenCleaner generates `Cleaned HPED(CHPED)` file, e.g. '\*.chped', which contains the converted HLA alleles and its log file, e.g. '\*.chped.log'.

You can find more information about the official HLA allele nomenclature defined by "IMGT-HLA" organization([go to IMGT-HLA Nomenclature](http://hla.alleles.org/nomenclature/naming.html)), which indicates that "Maximum 4-field Nomenclature" is currently used as the standard for HLA allele names.

<br>
When converting a lower resolution allele to a higher resolution (two-field to four-field, say), this module neither conduct a prediction nor make a guess. It just chooses a first option among possible candidates. For example,if  a HLA allele "A\*01:01" is given, then the possible candidates are

Allele | G_group | P_group
-------|---------|---------
A\*01:01:01:01|01:01:01G|01:01P
A\*01:01:01:02N|01:01:01G|0
A\*01:01:01:03|01:01:01G|01:01P
A\*01:01:01:04|01:01:01G|01:01P
A\*01:01:01:05|01:01:01G|01:01P
A\*01:01:01:06|01:01:01G|01:01P
...|...|...

The total number of possible candidates for "A\*01:01" is 96. But this module will just transform "A\*01:01" to "A\*01:01:01:01".

<br>

It is important to format the PED file correctly.
HLA PED file contains the HLA alleles of eight classical genes (FID,IID,pID,mID,SEX,PHENO,A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1).
HLA alleles must be in the following (alphabetical) order:
  HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1  

Note that each HLA allele spans two columns per individual (i.e. two chromosomes).

Put "0" for unknown HLA alleles.

As a result, an example row of HLA PED file will be

`
HG00319 HG00319 0       0       0       0       A*02:01:01G     A*03:01:01G     B*07:02:01G     B*56:01:01G     C*01:02:01G     C*03:04:01G     DPA1*01:03:01G  DPA1*01:03:01G  DPB1*04:01:01G  DPB1*04:01:01G  DQA1*01:01:01G  DQA1*01:02:01G  DQB1*05:01:01G  DQB1*06:02:01G  DRB1*01:01:01G  DRB1*15:01:01G`

Look at our example file example/g1k_subset.ped for an example formatting.

## (3) Usage examples

NomenCleaner in HLA-TAPAS has to be implemented in the directory of main project folder. (i.e. 'HLA-TAPAS/' where 'HLA-TAPAS.py' script is.)

```
$ cd ../ 
# Change your current directory to the HLA-TAPAS main project folder.
```

```
$ python -m NomenCleaner \
    --hped NomenCleaner/example/g1k_subset.ped \
    --out NomenCleaner/example/g1k_subset \
    --4field
```
