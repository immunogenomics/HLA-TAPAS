# MakeReference_v2

## (1) Introduction.

This is copy of project which is supposed to be sent to Y. Luo.

Please check 'docs/' directory for more detailed documentation("MakeReference", "NomenCleaner").


<br>
<br>


## (2) Usage example

MakeReference_v2 in HLA-TAPAS has to be implemented in the directory of main project folder. (i.e. 'HLA-TAPAS/' where 'HLA-TAPAS.py' script is.)

```
$ cd ../ 
# Change your current directory to the HLA-TAPAS main project folder.
```

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
