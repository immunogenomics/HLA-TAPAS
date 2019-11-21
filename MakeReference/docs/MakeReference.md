# MakeReference(2.0) 
------

Original Version Author: Sherman Jia (xiaomingjia@gmail.com)  
Recoded by W. Choi(wschoi.bhlab@gmail.com) and B. Han(buhm.han@gmail.com).


<br>
### [ Introduction ]

Original version of **MakeReference**([link](http://software.broadinstitute.org/mpg/snp2hla/makereference_manual.html)) was created by Sherman Jia. It generates a custom reference panel given HLA types and SNPs for HLA imputation. The panel works as indispensable part for **SNP2HLA**([link](http://software.broadinstitute.org/mpg/snp2hla/)).  

Original MakeReference works well, but it has 2 major shortcomings.
1. HLA allele name is restricted to "2-field(4-digits) Nomenclature System".
2. IMGT-HLA version is restricted to "IMGT/HLA Release: 3.7.0" version.  

The new version of **MakeReference(2.0)** is made to cope with these problems. It uses a ped file whose HLA alleles are in "Standard (maximum) 4-field Nomenclature System" instead of "2-field system"(ex. "HAPMAP_CEU_HLA.ped"). 
The advantage is that, if the input is in 3 or 4 fields from NGS typing, we don't have to lose accuracy by enforcing them to be 2 fields as in the old version of MakeReference. If the input is in 2 fields, the new version works exactly the same as the old version.



<br>
### [ Input and Output ]

Overall, Data files needed to generate custom reference panel are almost the same as those used in original **MakeReference**.

1. SNP dataset (ex. "HAPMAP_CEU.{bed,bim,fam}")
2. HLA types for SNP dataset individuals (ex. "HAPMAP_CEU_HLA.4field.ped")
3. HLA sequence Information (ex. "HLA_DICTIONARY_AA.hg18.imgt370.{txt,map}")

You can find above data files in 'data/MakeReference' directory which are prepared as an example data.

However, When it comes to HLA type data and HLA sequence information data, HLA allele name in these files became quite different. For example, content of ped file use in original **MakeReference** looks like this.

[//]: # (Original 2-field ped file content.)
FID|IID|pID|mID|SEX|PHENO|A_1|A_2|B_1|B_2|C_1|C_2
---|
1334|NA10847|NA12146|NA12239|2|0|2501|0301|0801|1801|0701|1203

\

DPA1_1|DPA1_2|DPB1_1|DPB1_2|DQA1_1|DQA_2|DQB1_1|DQB1_2|DRB1_1|DRB1_2
-|
0|0|0|0|0101|0101|0501|0501|0101|0101

But, the new MakeReference now uses a ped file like this.

[//]: # (New ped standard 4-field file content.)
FID|IID|pID|mID|SEX|PHENO|A_1|A_2|B_1|B_2
---|
1334|NA10847|NA12146|NA12239|2|0|A*25:01:01|A*03:01:01:01|B*08:01:01|B*18:01:01

\

C_1|C_2|DPA1_1|DPA1_2|DPB1_1|DPB1_2|DQA1_1|DQA_2
---|
C*07:01:01|C*12:03:01:01|0|0|0|0|DQA1*01:01:01|DQA1*01:01:01

\

DQB1_1|DQB1_2|DRB1_1|DRB1_2
------|
DQB1*05:01:01:01|DQB1*05:01:01:01|DRB1*01:01:01|DRB1*01:01:01

Also, HLA allele name in HLA sequence dictionary file will be prepared based on standard 4-field HLA allele name. For instance,

[//]: # (Before and After of content of HLA Sequence Dictionary file.)
(Before: "data/MakeReference_old/HLA_DICTIONARY_AA.txt")

Allele|Seq
------|---
A:01:01|MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGR...

(After: "data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt")

Allele|Seq
------|---
A*01:01:01:01|MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGR...



<br>
> cf) A ped file in the old style can be converted to the new 'standard 4-field' by "NomenCleaner.py" module that we provide.




<br>
### [ Software Requirements ]

#### (External Software)

1.  beagle.jar(v.3.0.4)
2.  beagle2vcf.jar
3.  linkage2beagle.jar
4.  plink(v.1.07)

These software must be prepared in "dependency/" directory.

#### (Python)

**MakeReference(2.0)** works primarily based on 
> "Python 3.6.x :: Anaconda custom (x86_64)" (OS X)  

> "Python 3.6.x :: Anaconda, Inc." (Linux(CentOS))

I stronlgy recommend for users installing Anaconda with Python 3.6 version.

#### (Perl)

> Perl v5.10.x to v5.18.x recommended.

#### (Java)

> Java 1.8.x to 10.0.1 version recommended.


<br>
### [ Example ]


<br>
#### (Case 1): MakeReference(2.0) / imgt3.7.0 / hg18  

``` console

$ python3 MakeReference.py \
          -i ./data/MakeReference/HAPMAP_CEU \
          -ped data/MakeReference/HAPMAP_CEU_HLA.4field.ped \
          -hg 18 \
          -o ./NEWVERSION_imgt370_hg18/HAPMAP_CEU \
          -dict-AA data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt \
          -dict-AA-map ./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map \
          -dict-SNPS ./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt \
          -dict-SNPS-map ./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map
        
```

Users can obtain the same result of "(Case 1)" with alleles in standard 4-field HLA allele name. The only difference is that **MakeReference(2.0)** uses .ped file with 4-field and HLA sequence information dictionaries should be given through those 4 arguments("-dict-AA", ...).


<br>
#### (Case 2): MakeReference(2.0) / imgt3.32.0 / hg19

``` console

$ python3 MakeReference.py \
          -i data/MakeReference/HAPMAP_CEU \
          -ped data/MakeReference/HAPMAP_CEU_HLA.4field.ped \
          -hg 19 \
          -o NEWVERSION_imgt3320_hg19 \
          -dict-AA ./data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.txt \
          -dict-AA-map ./data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.map \
          -dict-SNPS ./data/MakeReference/HLA_DICTIONARY_SNPS.hg19.imgt3320.txt \
          -dict-SNPS-map ./data/MakeReference/HLA_DICTIONARY_SNPS.hg19.imgt3320.map

```



#### (Case 3): Old version usage. 
**MakeReference(2.0)** can work exactly the same as original **MakeReference** given "--previous-version" argument.

``` console

$ python3 MakeReference.py \
          --previous-version \
          -i ./data/MakeReference/HAPMAP_CEU \
          -ped data/MakeReference_old/HAPMAP_CEU_HLA.ped \
          -hg 18 \
          -o ./NEWVERSION_imgt370_hg18/HAMAP_CEU
          
```



> (Note) Original **MakeReference** only supports hg 18. So, if you give "-hg 19" or "-hg 38" argument with "--previous-version", then it will be overriden to "-hg 18".
