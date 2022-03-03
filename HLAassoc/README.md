# HLAassoc

## (1) Introduction
This tapa provides a HLA-focused statistical association analysis. Depends on the disease of focus, we can either apply a linear (quantitative trait) or logisitc (binary trait) for single variant association test. We can also perform a conditional analysis (`--condition` or `--condition-list`) for fine-mapping studies.

Due to the high LD structure within the extended MHC region, a **haplotype** based (`OMNIBUS`) approach is often implemented in the HLA fine-mapping studies. For each amino acid position, we perform a multiallelic association between trait of interest and a haplotype matrix (of the amino acid positions). To get an omnibus P-value for each position, we estimated the effect of each amino acid by assessing the significance of the improvement in fit by calculating the in-model fit, compared to a null model following an F-distribution with m-1 degrees of freedom, where m is the number of residues in an amino acid position. This is implemented using an ANOVA test in R. The most frequent haplotype was excluded from a haplotype matrix as a reference haplotype for association. 

<br>

## (2) Input and output
You would need the imputed `vcf` file and a phenotype file (`--pheno`) to start your association test. You can alsp provide other covariates using the `--covar` flag.


## (3) Usage examples

HLAassoc in HLA-TAPAS has to be implemented in the directory of main project folder. (i.e. 'HLA-TAPAS/' where 'HLA-TAPAS.py' script is.)

```
$ cd ../ 
# Change your current directory to the HLA-TAPAS main project folder.
```

<br>

### (a) Logistic Regression

```
$ python -m HLAassoc LOGISTIC \
    --vcf HLAassoc/example/IMPUTED.1958BC.bgl.phased.vcf.gz \
    --out Myassoc/IMPUTED.1958BC \
    --pheno HLAassoc/example/1958BC.phe \
    --pheno-name p1
```

If user wants to make HLA marker have un-transformed allele name, say if your reference panel is only at the G-group resolution, you might want to transform the 4-field alleles back to G-group alleles for association studies. Users can get this by passing '--hped' and '--chped' arguments.

```
$ python -m HLAassoc LOGISTIC \
    --vcf HLAassoc/example/IMPUTED.1958BC.bgl.phased.vcf.gz \
    --out Myassoc/IMPUTED.1958BC.rev_map \
    --pheno HLassoc/example/1958BC.phe \
    --pheno-name p1 \
    --hped HLAassoc/example/1958BC.Ggroup.hped \
    --chped HLAassoc/example/1958BC.imgt3320.4field.chped
```
<br>

### (b) Linear Regression 

Similar to the logisitc regression, here's how linear regression works:

```
$ python -m HLAassoc LINEAR \
    --vcf HLAassoc/example/IMPUTED.1958BC.bgl.phased.vcf.gz \
    --out Myassoc/IMPUTED.1958BC \
    --pheno HLAassoc/example/1958BC.phe \
    --pheno-name p2
```
<br>

### (c) Omnibus Test

An omnibus test is performed to determine the associate between each **amino acid position** and phenotype of interest, accounting for all _m_ amnio acid residue changes occuring at that position. This test is equivalent to testing a single multi-allelic locus for association with _m_ alleles. 


Assuming you have an imputed vcf file (e.g., from the Michigan imputation server), next minimal 3 files can be passed to perform an omnibus test.

1. (Imputed) VCF file generated from SNP2HLA. (*.vcf.gz)
2. PLINK bim file of reference panel. (*.bim)
3. PLINK fam file of target samples. (*.fam)

PLINK fam file must Phenotype information. All samples will be considered to be originated from the same population.

```
$ python -m HLAassoc OMNIBUS \
    --vcf HLAassoc/example/OMNIBUS/Case+Control+1000G_EUR_REF.IMPUTED.bgl.phased.vcf.gz \
    --bim HLAassoc/example/OMNIBUS/Case+Control+1000G_EUR_REF.REF.bglv4.bim \
    --fam HLAassoc/example/OMNIBUS/Case+Control.300+300.chr6.hg18.fam \
    --out Myassoc/Case+Control+1000G_EUR_REF.OMNIBUS \
    --aa-only \
    --maf-threshold 0
```

<!-- <br>

```
$ python -m HLAassoc OMNIBUS \
    --vcf HLAassoc/example/OMNIBUS/1958BC+HM_CEU_REF.IMPUTED.bgl.phased.vcf.gz \
    --bim HLAassoc/example/OMNIBUS/1958BC+HM_CEU_REF.REF.bglv4.bim \
    --fam HLAassoc/example/OMNIBUS/1958BC.fam \
    --pheno HLAassoc/example/OMNIBUS/1958BC.phe \
    --out Myassoc/1958BC+HM_CEU_REF.minimal_input.OMNIBUS \
    --aa-only \
    --maf-threshold 0
``` -->

Omnibus Test requires next 5 files, each of which contains necessary information to perform an omnibus test.

1. (Imputed) Phased BEAGLE file. (*.bgl.phased)
2. PLINK bim file of **reference** panel. (*.bim)
3. PLINK fam file of target samples. (*.fam)
4. A covariate file e.g. genetic principal components, sex, age, etc. (*.covs; with header: FID, IID, sex, age, PC1, PC2, etc..)
5. Phenotype file of target samples. (*.pheno;Case: 1/ Control: 0; **no** header. Three columns: FID, IID, Pheno)

Assuming 1 - 5 files have the same file prefix and those specified file extensions, those 5 files can be passed to the Omnibus Test **all at once** with the '--file' argument.

<!-- ```
# Under construction. (2020.05.11. WS.)
$ Rscript run_omnibus_test.R --file cohort --pop pop  \
	--out output --aa-only --omnibus \
	--remove-samples-by-haplo \
     	--remove-samples-aa-pattern AA_B \
	--min-haplo-count 10 --maf-threshold 0 
``` -->

```
$ python -m HLAassoc OMNIBUS \
    --file HLAassoc/example/OMNIBUS/Case+Control+1000G_EUR_REF.IMPUTED.chr6.hg18.100+100 \
    --pop HLAassoc/example/OMNIBUS/Case+Control+1000G_EUR_REF.IMPUTED.chr6.hg18.100+100.pop \
    --out Myassoc/Case+Control+1000G_EUR_REF.OMNIBUS \
    --aa-only \
    --maf-threshold 0
```

The above command is equivalent to the next command.

```
$ python -m HLAassoc OMNIBUS \
    --phased HLAassoc/example/OMNIBUS/Case+Control+1000G_EUR_REF.IMPUTED.chr6.hg18.100+100.bgl.phased \
    --fam HLAassoc/example/OMNIBUS/Case+Control+1000G_EUR_REF.IMPUTED.chr6.hg18.100+100.fam \
    --covars HLAassoc/example/OMNIBUS/Case+Control+1000G_EUR_REF.IMPUTED.chr6.hg18.100+100.pcs \
    --pheno HLAassoc/example/OMNIBUS/Case+Control+1000G_EUR_REF.IMPUTED.chr6.hg18.100+100.pheno \
    --out Myassoc/Case+Control+1000G_EUR_REF.OMNIBUS \
    --aa-only \
    --maf-threshold 0
```

<br>

