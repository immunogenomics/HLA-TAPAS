# Input
1. SNP dataset (.bed/bim/fam PLINK format)
2. Reference dataset (.bgl.phased.vcf.gz/.markers Beagle)

# Running command
	`python3 SNP2HLA.py -i DATA(.bed/.bim/.fam) -o OUTPUT -rf REFERENCE(.bgl.phased/.markers) --nthreads NTHREAD --java-mem=MEMORY[GB]`

# Example
Run SNP2HLA with sample data provided (10 samples) using the following command:

	python3 SNP2HLA.py -i data/test -o test -rf data/T1DGCb37.bglv4 --nthreads 16 --java-mem=8g
