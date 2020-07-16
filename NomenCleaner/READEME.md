# NomenCleaner

## (1) Introduction

`NomenCleaner` is a module to help researchers to solve challenges related to HLA nomenclature. Two major challenges would be (1) Conversion between nomenclatures, and (2) Field checking.



1. Conversion

NomenCleaner converts a given set of HLA alleles which is in arbitrary nomenclature to the ones in the updated nomenclature. **The conversion from the updated to old nomenclature is NOT available**.

<br>

2. Field checking

**Often, the field separator is removed though given HLA allele conforms to the updated nomenclature**. This barely causes a problem in most cases, however, there is definitely some case that may confuse researchers. For example, if the HLA allele **DPB1\*101101** is given, then some researchers may be perplexed as to which one is right among <U>DBP1\*10:11:01</U>, <U>DPB1\*101:101</U> and <U>DPB1\*1011:01</U>. NomenCleaner searches the `HAT(HLA Allele Table)` file, which contains whole HLA allele name information of the IMGT database of a specific version, and determines that **DPB1*1011:01 is the only valid solution**. (cf. DBP1*1011:01 is in the IMGT database of the version 3.37.0.)


<br>
<br>


## (2) Usage examples

NomenCleaner basically takes `HPED` file(e.g. '\*.hped') and `HAT(HLA Allele Table)` file (e.g. '\*.hat') as input, then performs its jobs as required by the user. Finally the NomenCleaner generates `Cleaned HPED(CHPED)` file, e.g. '\*.chped', which contains the converted HLA alleles and its log file, e.g. '\*.chped.log'.

NomenCleaner in HLA-TAPAS has to be implemented in the directory of main project folder. (i.e. 'HLA-TAPAS/' where 'HLA-TAPAS.py' script is.)

```
$ cd ../ 
# Change your current directory to the HLA-TAPAS main project folder.
```

```
$ python -m NomenCleaner \
    --hped NomenCleaner/example/wtccc_filtered_58C_RA.hatk.300+300.hped \
    --out RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.Ggroup \
    --Ggroup
```