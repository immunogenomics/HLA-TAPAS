# HLA-focused Manhattan Plot

## (1) Introduction
This tapa provides basic functions for a HLA-focused Manhattan plot. This module is adapted from [HATK](https://github.com/WansonChoi/HATK). The script searches for all association results in a directory (e.g. *.assoc or *.logistic.assoc) and plot all association results in a single figure (see below)

![Conditional Manhattan Plot](../../manuscript/figs/SF18_all_conditional.png)

<br>
<br>

## (2) Usage example

The module 'Manhattan' in HLA-TAPAS has to be implemented in the directory of main project folder. (i.e. 'HLA-TAPAS/' where 'HLA-TAPAS.py' script is.) and the user needs to specify the association folder.

```
$ cd ../ 
# Change your current directory to the HLA-TAPAS main project folder.
```

```
$ python -m HLAManhattan \
    --assoc-result HLAManhattan/example/1958BC+HM_CEU_REF.IMPUTED.assoc.logistic \
    --hg 18 \
    --out MyManhattan/1958BC+HM_CEU_REF.IMPUTED
```
