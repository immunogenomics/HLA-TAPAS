# SNP2HLA

## (1) Introduction
SNP2HLA module (tapa) performs HLA imputation (amino acids, *HLA* alleles and SNPs) using genotyping data.

## (2) Input & output

You would need two things for HLA imputation:

1. A GWAS dataset contains the MHC region (.bed/bim/fam PLINK format)
**N.B.: due to the BEAGLE update, unlike SNP2HLA-v1, the correct genomic build (hg19/Grch38) is now important for imputation**
2. Reference dataset (.bgl.phased.vcf.gz/.markers/.FRQ.frq/.bim; output from the [MakeReference](../MakeReference) module)

After imputaton, SNP2HLA will output a single vcf file (\*.bgl.phased.vcf.gz).


## (2) Usage example

SNP2HLA in HLA-TAPAS has to be implemented in the directory of main project folder. (i.e. 'HLA-TAPAS/' where 'HLA-TAPAS.py' script is.)

```
$ cd ../ 
# Change your current directory to the HLA-TAPAS main project folder.
```
```
$ python3 -m SNP2HLA \
	--target SNP2HLA/example/1958BC \
	--out SNP2HLA/example/IMPUTED.1958BC \
	--reference resources/1000G.bglv4 \
	--nthreads 2 \
	--mem 4g
```
