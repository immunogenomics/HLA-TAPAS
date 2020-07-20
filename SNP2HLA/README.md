# SNP2HLA

## (1) Introduction)
1. SNP dataset (.bed/bim/fam PLINK format)
2. Reference dataset (.bgl.phased.vcf.gz/.markers; Beagle4.1)

## (2) Usage example

SNP2HLA in HLA-TAPAS has to be implemented in the directory of main project folder. (i.e. 'HLA-TAPAS/' where 'HLA-TAPAS.py' script is.)

```
$ cd ../ 
# Change your current directory to the HLA-TAPAS main project folder.
```
```
$ python3 -m SNP2HLA \
	--target SNP2HLA/example/1958BC \
	--out MySNP2HLA/IMPUTED.1958BC \
	--reference SNP2HLA/example/T1DGCb37.bglv4 \
	--nthreads 2 \
	--mem 4g
```
