#!/bin/csh -f

######################################################################################################
#
# SNP2HLA_new.csh
#
# Author: Sherman Jia (xiaomingjia@gmail.com)
#         + Small modifications by Buhm Han (buhmhan@broadinstitute.org): 8/7/12
#         + Extensive modifications by Phil Stuart (pstuart@umich.edu) on 4/19/16 to allow use of Beagle 4.1:
#           verified to work with 22Feb16, 15Apr16, and 03May16 versions of Beagle.
#         + Small modifications by Yang Luo (yangluo@broadinstitute.org): 09/30/16: verfiied working with Bealge 4.1 27Jun16.
# DESCRIPTION: This script runs imputation of HLA amino acids and classical alleles using SNP data.
#
# INPUTS:
# 1. Plink dataset (*.bed/bim/fam)
# 2. Reference dataset (*.bgl.phased, *.markers in beagle 3.0.4 format; *.fam/.bim/.FRQ.frq in PLINK format)
#
# DEPENDENCIES: (download and place in the same folder as this script)
# 1. PLINK (1.9)  (Will not work with older Plink 1.07)
# 2. Beagle (4.1) (Need to rename java executable as beagle.jar)
# 3. merge_tables.pl (Perl script to merge files indexed by a specific column)
# 4. vcf2gprobs.jar (Beagle utility for generating a Beagle v3 genotypes probability file from a Beagle 4.1 vcf file with GT field data)
# 5. beagle2vcf and vcf2phased (Utilities by Phil Stuart for vcf <-> Beagle v3 format)
# 6. ParseDosage.csh (Converts Beagle posterior probabilities [.gprobs] to dosages in PLINK format [.dosage])
# 7. revert_alleles (Utility by Phil Stuart to revert hacked non-ACTG alleles in *.vcf, *.gprobs and *.dosage output files (necessitated by Beagle 4.1) back to allele identities in reference panel)
# 8. If genetic_map_file argument is specified, PLINK format genetic map on cM scale (plink.chr6.GRCh36.map, downloaded from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)
#
# USAGE: ./SNP2HLA_new.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers/.fam/.bim/.bed/.FRQ.frq) OUTPUT plink max_memory[gb] nthreads niterations genetic_map_file
#
######################################################################################################

if ($#argv < 4) then
    echo "USAGE: ./SNP2HLA_new.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers/.fam/.bim/.bed/.FRQ.frq) OUTPUT plink {optional: java_max_memory[gb] number_of_threads number_of_iterations plink_format_genetic_map_file}"; exit 1
endif

set SCRIPTPATH=`dirname $0`

set MERGE=$SCRIPTPATH/merge_tables.pl
set PARSEDOSAGE=$SCRIPTPATH/ParseDosage.csh

# CHECK FOR DEPENDENCIES
if (! -e `which $4`) then
    echo "Please install version 1.9 of PLINK (https://www.cog-genomics.org/plink2) and point to the plink run file.";
    echo "tcsh: use plink"
    echo "bash: use ./plink"
    exit 1
else if (! -e $SCRIPTPATH/beagle.jar) then
    echo "Please install Beagle 4.1 (https://faculty.washington.edu/browning/beagle/beagle.html#download) and rename the run file as beagle.jar and copy it into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/vcf2gprobs.jar) then
    echo "Please copy Beagle utility file vcf2gprobs.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html#vcf2beagle) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/beagle2vcf) then
    echo "Please copy beagle2vcf (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/vcf2phased) then # We use beagle2linkage (Buhm, 8/13/12)
    echo "Please copy vcf2phased (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $MERGE) then
    echo "Please copy merge_tables.pl (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $PARSEDOSAGE) then
    echo "Please copy ParseDosage.csh (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/revert_alleles) then
    echo "Please copy revert_alleles (included with this package) into $SCRIPTPATH/"; exit 1
endif

# INPUTS
set INPUT=$1
set REFERENCE=$2
set OUTPUT=$3
set PLINK=$4

if ($#argv >= 5) then
    set MEM=$5
else
    set MEM=2 # Default java memory 2Gb
endif

if ($#argv >= 6) then
    set THREAD=$6
else
    set THREAD=1
endif

if ($#argv >= 7) then
    set ITER=$7
else
    set ITER=5
endif

if ($#argv >= 8) then
    set MAP=$8
endif

#CHECK FOR PRESENCE OF INPUT FILES IN DIRECTORY

if (! -f $INPUT.bed) then
  echo "Input file $INPUT.bed not found"; exit 1
else if (! -f $INPUT.bim) then
  echo "Input file $INPUT.bim not found"; exit 1
else if (! -f $INPUT.fam) then
  echo "Input file $INPUT.fam not found"; exit 1
else if (! -f $INPUT.bim) then
  echo "Input file $INPUT.bed not found"; exit 1
else if (! -f $REFERENCE.bim) then
  echo "Input file $REFERENCE.bim not found"; exit 1
else if (! -f $REFERENCE.markers) then
  echo "Input file $REFERENCE.markers not found"; exit 1
else if (! -f $REFERENCE.fam) then
  echo "Input file $REFERENCE.fam not found"; exit 1
else if (! -f $REFERENCE.bgl.phased) then
  echo "Input file $REFERENCE.bgl.phased not found"; exit 1
else if (! -f $REFERENCE.FRQ.frq) then
  echo "Input file $REFERENCE.FRQ.frq not found"; exit 1
endif

set JAVATMP=$OUTPUT.javatmpdir
mkdir -p $JAVATMP
alias plink '$PLINK --noweb --silent --allow-no-sex'
alias beagle 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\g -jar $SCRIPTPATH/beagle.jar'
alias vcf2gprobs 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\g -jar $SCRIPTPATH/vcf2gprobs.jar'
alias beagle2vcf '$SCRIPTPATH/beagle2vcf -chr 6 -missing 0 -recode 1 -gtype_opt 1'
alias vcf2phased '$SCRIPTPATH/vcf2phased'
alias revert_alleles '$SCRIPTPATH/revert_alleles'

# Functions to run
set CONVERT_IN  = 0 #no need to convert again using the new beagle version
set EXTRACT_MHC = 1
set FLIP        = 1
set IMPUTE      = 1
set CONVERT_OUT = 0 #convert back to beagle v3 formates
set CLEANUP     = 1

# SET PARAMETERS
set TOLERATED_DIFF = .15
set i = 1

echo ""
echo "SNP2HLA: Performing HLA imputation for dataset $INPUT";
echo "- Java memory = "$MEM"Gb"
echo "- Beagle number of threads = "$THREAD" threads"
echo "- Beagle number of phasing iterations = "$ITER" iterations"
if ($#argv >= 8) then
    echo "- Beagle genetic map file = "$MAP""
else
    echo "- Beagle genetic map file = none specified"
endif

# Convert reference panel to vcf format required by Beagle 4.1

if ($CONVERT_IN) then
    echo ""
    echo "[$i] Converting reference panel to vcf format required by Beagle 4.1."; @ i++
    ln -s $REFERENCE.markers temp.markers
    ln -s $REFERENCE.fam temp.fam
    ln -s $REFERENCE.bgl.phased temp.bgl.phased
    beagle2vcf -fnmarker temp.markers -fnfam temp.fam -fngtype temp.bgl.phased -fnout $REFERENCE.vcf
    rm temp.markers
    rm temp.fam
    rm temp.bgl.phased
endif

set MHC=$OUTPUT.MHC

if ($EXTRACT_MHC) then
    echo "[$i] Extracting SNPs from the MHC."; @ i++
    plink --bfile $INPUT --chr 6 --from-mb 29 --to-mb 34 --maf 0.01 --make-bed --out $OUTPUT.MHC #be less restrict about maf 1%
endif

if ($FLIP) then
    echo ""
    echo "[$i] Performing SNP quality control."; @ i++

    # Identifying non-A/T non-C/G SNPs to flip
    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.bim >> $OUTPUT.tmp1
    echo "SNP 	POSR	A1R	A2R" > $OUTPUT.tmp2
    cut -f2,4- $REFERENCE.bim >> $OUTPUT.tmp2
    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP |  grep -v -w NA > $OUTPUT.SNPS.alleles

    awk '{if ($3 != $6 && $3 != $7){print $1}}' $OUTPUT.SNPS.alleles > $OUTPUT.SNPS.toflip1
    plink --bfile $MHC --flip $OUTPUT.SNPS.toflip1 --make-bed --out $MHC.FLP

    # Calculating allele frequencies
    plink --bfile $MHC.FLP --freq --out $MHC.FLP.FRQ
    sed 's/A1/A1I/g' $MHC.FLP.FRQ.frq | sed 's/A2/A2I/g' | sed 's/MAF/MAF_I/g' > $OUTPUT.tmp

    mv $OUTPUT.tmp $MHC.FLP.FRQ
    $MERGE $REFERENCE.FRQ.frq $MHC.FLP.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.frq
    sed 's/ /\t/g' $OUTPUT.SNPS.frq | awk '{if ($3 != $8){print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $8 "\t" 1-$10 "\t*"}else{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 "\t."}}' > $OUTPUT.SNPS.frq.parsed

    # Finding A/T and C/G SNPs
    awk '{if (($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C")){if ($4 > $7){diff=$4 - $7; if ($4 > 1-$7){corrected=$4-(1-$7)}else{corrected=(1-$7)-$4}}else{diff=$7-$4;if($7 > (1-$4)){corrected=$7-(1-$4)}else{corrected=(1-$4)-$7}};print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" diff "\t" corrected}}' $OUTPUT.SNPS.frq.parsed > $OUTPUT.SNPS.ATCG.frq

    # Identifying A/T and C/G SNPs to flip or remove
    awk '{if ($10 < $9 && $10 < .15){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toflip2
    awk '{if ($4 > 0.4){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toremove

    # Identifying non A/T and non C/G SNPs to remove
    awk '{if (!(($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C"))){if ($4 > $7){diff=$4 - $7;}else{diff=$7-$4}; if (diff > '$TOLERATED_DIFF'){print $1}}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    awk '{if (($2 != "A" && $2 != "C" && $2 != "G" && $2 != "T") || ($3 != "A" && $3 != "C" && $3 != "G" && $3 != "T")){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    awk '{if (($2 == $5 && $3 != $6) || ($3 == $6 && $2 != $5)){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove

    # Making QCd SNP file
    plink --bfile $MHC.FLP --geno 0.2 --exclude $OUTPUT.SNPS.toremove --flip $OUTPUT.SNPS.toflip2 --make-bed --out $MHC.QC
    plink --bfile $MHC.QC --freq --out $MHC.QC.FRQ
    sed 's/A1/A1I/g' $MHC.QC.FRQ.frq | sed 's/A2/A2I/g' | sed 's/MAF/MAF_I/g' > $OUTPUT.tmp
    mv $OUTPUT.tmp $MHC.QC.FRQ.frq
    $MERGE $REFERENCE.FRQ.frq $MHC.QC.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.QC.frq

    cut -f2 $OUTPUT.SNPS.QC.frq | awk '{if (NR > 1){print $1}}' > $OUTPUT.SNPS.toinclude

    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.QC.bim >> $OUTPUT.tmp1

    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP | awk '{if (NR > 1){if ($5 != "NA"){pos=$5}else{pos=$2}; print "6\t" $1 "\t0\t" pos "\t" $3 "\t" $4}}' > $MHC.QC.bim

    # Extracting SNPs and recoding QC'd file as vcf
    plink --bfile $MHC.QC --extract $OUTPUT.SNPS.toinclude --make-bed --out $MHC.QC.reorder
	plink --bfile $MHC.QC.reorder --recode vcf-iid --a1-allele $REFERENCE.markers 4 1 --out $MHC.QC

    # Remove temporary files
    rm $OUTPUT.tmp1 $OUTPUT.tmp2
    rm $MHC.FLP*
    rm $OUTPUT.SNPS.*
endif

if ($IMPUTE) then
    echo ""
    echo "[$i] Performing HLA imputation (see $OUTPUT.bgl.log for progress)."; @ i++

	if ($#argv >= 8) then
        beagle ref=$REFERENCE.bgl.vcf.gz gt=$MHC.QC.vcf impute=true gprobs=true nthreads=$THREAD chrom=6 niterations=$ITER lowmem=true out=$OUTPUT.bgl map=$MAP
	else
        beagle ref=$REFERENCE.bgl.vcf.gz gt=$MHC.QC.vcf impute=true gprobs=true nthreads=$THREAD chrom=6 niterations=$ITER lowmem=true out=$OUTPUT.bgl
	endif
endif

if ($CONVERT_OUT) then
    echo ""
    echo "[$i] Converting posterior probabilities to PLINK dosage format."; @ i++

    # Converting .gprobs to .dosage format
	gunzip -d $OUTPUT.bgl.vcf.gz
	cat $OUTPUT.bgl.vcf | vcf2gprobs > $OUTPUT.bgl.gprobs
	cp $OUTPUT.bgl.gprobs temp.gprobs
	sed -i '1d' temp.gprobs
	$PARSEDOSAGE temp.gprobs > $OUTPUT.dosage
	rm temp.gprobs

  echo ""
	echo "[$i] Converting vcf output to Beagle v3 phased haplotype and imputation quality files."; @ i++
  ln -s $REFERENCE.markers temp.markers
  ln -s $INPUT.fam temp.fam
  #set marker=`basename $REFERENCE.markers`
	vcf2phased -gprobs_opt 1 -fnvcf $OUTPUT.bgl.vcf -fnfam temp.fam -fnmarker temp.markers -fnphased $OUTPUT.bgl.phased -fnr2 $OUTPUT.bgl.r2
  rm temp.fam
  rm temp.markers

    echo ""
	echo "[$i] Reverting hacked allele identities in *.vcf, *.gprobs and *.dosage output files back to those used in reference panel."; @ i++

  ln -s $REFERENCE.bim temp.bim
  #set bim=`basename $REFERENCE.bim`
	revert_alleles -script_opt 1 -fnvcf $OUTPUT.bgl.vcf -fngprobs $OUTPUT.bgl.gprobs -fndosage $OUTPUT.dosage -fnref temp.bim
  rm temp.bim

endif

if ($CLEANUP) then
    rm $OUTPUT.MHC.*
    rm $OUTPUT.tmp*
    rm $OUTPUT.IMPUTED.*.bgl.phased.phased
    rm -r $JAVATMP
    rm -f plink.log
    rm *.nosex
    echo "DONE!"
    echo ""
endif
