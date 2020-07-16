# HLA-TAPAS
**HLA-TAPAS** (HLA-Typing At Protein for Association Studies) is an HLA-focused pipeline that can handle HLA reference panel construction (*MakeReference*), HLA imputation (*SNP2HLA*), and HLA association (*HLAassoc*). 
It is an updated version of the [SNP2HLA](http://software.broadinstitute.org/mpg/snp2hla/). 
Briefly, major updates include 

(1) using PLINK-v1.9 instead of v1.07; 

(2) using BEAGLE v4.1 (URLs) instead of v3 for phasing and imputation; and 

(3) including custom R scripts for performing association and fine-mapping analysis in multiple ancestries. 

## Requirments & Dependencies

Linux or OS_X environment is required for HLA-TAPAS. HLA-TAPAS currently doesn't support Windows.

Python system requires next settings.
- python=3.7.x
- pandas=1.0.3

R statistical programming language requires next settings.
- R=3.6.x
- argparse
- stringr
- purrr
- dplyr
- multidplyr
- tidyr
- data.table
- parallel
- rcompanion

Also, Next external software have to be prepared in 'dependency/' folder.
- PLINK v1.9b (https://www.cog-genomics.org/plink2)
- BEAGLE v4.1 (https://faculty.washington.edu/browning/beagle/b4_1.html#download)
- beagle2vcf.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html)
- linkage2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html)
- vcf2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html)

Related to BEAGLE v4.1, after downloading and preparing it in 'dependency/' folder, **RENAME IT TO 'beagle.jar'**.

(Sorry Users! You need to download them yourselves due to copyright issue.)


<br>
<br>


## Usage example

Main usage example of HLA-TAPAS is as below.

```
$ python HLA-TAPAS.py \
    --target example/wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18 \
    --reference example/1000G.EUR.chr6.hg18.28mb-35mb \
    --hped-Ggroup example/1000G.EUR.Ggroup.hped \
    --pheno example/wtccc_filtered_58C_RA.hatk.300+300.phe \
    --hg 18 \
    --out MyHLA-TAPAS/WTCCC+1000G_EUR_REF \
#   --mem 4g \
#   --nthreads 4
```
Last two arguments are for Java Heap memory size and the number of threads to be used in Beagle. Users can specify those values on their own.

Each main module of HLA-TAPAS can be implemented separately.

(1) NomenCleaner
```
$ python -m NomenCleaner \
    --hped NomenCleaner/example/wtccc_filtered_58C_RA.hatk.300+300.hped \
    --out RESULT_EXAMPLE_wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.Ggroup \
    --Ggroup
```

(2) MakeReference
```
$ python -m MakeReference \
    --variants MakeReference/example/HAPMAP_CEU \
    --chped MakeReference/example/HAPMAP_CEU_HLA.imgt3320.4field.chped \
    --hg 19 \
    --out MyRef/HAPMAP_CEU.REF.bglv4 \
    --dict-AA MakeReference/data/hg19/HLA_DICTIONARY_AA.hg19.imgt3320 \
    --dict-SNPS MakeReference/data/hg19/HLA_DICTIONARY_SNPS.hg19.imgt3320 \
    --phasing
```

(3) SNP2HLA
```
$ python -m SNP2HLA \
    --target SNP2HLA/example/1958BC \
    --reference SNP2HLA/example/T1DGCb37.bglv4 \
    --out MySNP2HLA/IMPUTED.1958BC \
    --nthreads 2 \
    --mem 4g
```

(4) HLA-assoc
```
$ python -m HLA_assoc \
    --target HLA_assoc/example/IMPUTED.1958BC.bgl.phased.vcf.gz \
    --out MyHLAassoc/IMPUTED.1958BC \
    --pheno example/1958BC.phe \
    --pheno-name p1
```

(5) Manhattan
```
$ python -m Manhattan \
    --assoc-result Manhattan/example/1958BC+HM_CEU_REF.IMPUTED.assoc.logistic \
    --hg 18 \
    --out MyManhattan/1958BC+HM_CEU_REF.IMPUTED
```

<br>
<br>


<!-- 
## Development Log

(2020.04.30.) 

[MakeReference_v2]
    
- refined dependency checkings.
- introduced two tricks, (1) ATtrick and (2) redefineBP, to make Beagle framework work with binary markers.
- File conversion from PLINK to VCF format in 'PREPARE' code block.
- Phasing with Beagle 4.1 is now available.


[SNP2HLA]
- modified the way to convert target PLINK SNP data to VCF file.
- modified bash command for imputation and its execution way(os.system() -> subprocess.run()).

[NomenCleaner]
- disjointly integrated NomenCleaner as a package.
- Upgraded parts have been applied.


(2020.05.01)

[HLA_assoc]
- introduced reverse-mapping module which Yang requested in the recent mail(2020.04.26 00:22).
- introduced PLINK logistic regression for association test.

(2020.05.02.)

[HLA-TAPAS]
- All (1) NomenCleaner, (2) MakeReference(v2, Beagle4.1), (3) SNP2HLA(Python version, Beagle4.1) and (4) HLA_assoc have been integrated in 'HLA-TAPAS.py' main pipeline.
- succeeded in the 1st whole implementations of those functions with 'example/' data.
 -->
