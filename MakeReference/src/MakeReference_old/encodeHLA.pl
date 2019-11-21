#!/usr/bin/perl

###########################################################################################
#
# Sherman Jia, 2012
# encodeHLA.pl
#
# This script encodes HLA alleles into bi-allelic markers (for imputation and PLINK analysis)
# The input ped file should contain: FID,IID,pID,mID,SEX,PHENO,
#                                    2 each of: HLA-A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1
#
###########################################################################################

if (scalar(@ARGV) !=2){ die "usage: ./encodeHLA.pl HLApedFile outputmapname > newped\n";}

$pedname=shift;
$mapname=shift;

my $nsamples = 0;
my %alleles = ();
my %genepos = ();

# $genepos{"HLA_A"} =    30019970;
# $genepos{"HLA_C"} =    31346171;
# $genepos{"HLA_B"} =    31431272;
# $genepos{"HLA_DRB1"} = 32660042;
# $genepos{"HLA_DQA1"} = 32716284;
# $genepos{"HLA_DQB1"} = 32739039;
# $genepos{"HLA_DPA1"} = 33145064;
# $genepos{"HLA_DPB1"} = 33157346;

# new genomic positions 
$genepos{"HLA_A"} =    30019970 - 107979;
$genepos{"HLA_C"} =    31346171 - 107979;
$genepos{"HLA_B"} =    31431272 - 107979;
$genepos{"HLA_DRB1"} = 32660042 - 107979;
$genepos{"HLA_DQA1"} = 32716284 - 107979;
$genepos{"HLA_DQB1"} = 32739039 - 107979;
$genepos{"HLA_DPA1"} = 33145064 - 107979;
$genepos{"HLA_DPB1"} = 33157346 - 107979;

open(FILE,"$pedname") || die "can't open $pedname\n";
while (<FILE>) {
    #Adds each HLA allele to the hash table
    chomp;
    s/^\s+//;
    @fields = split(/\s+/);

    $alleles{ "HLA_A_".$fields[6] } = $genepos{"HLA_A"};     $alleles{ "HLA_A_".$fields[7] } = $genepos{"HLA_A"};
    $alleles{ "HLA_B_".$fields[8] } = $genepos{"HLA_B"};     $alleles{ "HLA_B_".$fields[9] } = $genepos{"HLA_B"};
    $alleles{ "HLA_C_".$fields[10] } = $genepos{"HLA_C"};    $alleles{ "HLA_C_".$fields[11] } = $genepos{"HLA_C"};
    $alleles{ "HLA_DPA1_".$fields[12] } = $genepos{"HLA_DPA1"}; $alleles{ "HLA_DPA1_".$fields[13] } = $genepos{"HLA_DPA1"};
    $alleles{ "HLA_DPB1_".$fields[14] } = $genepos{"HLA_DPB1"}; $alleles{ "HLA_DPB1_".$fields[15] } = $genepos{"HLA_DPB1"};
    $alleles{ "HLA_DQA1_".$fields[16] } = $genepos{"HLA_DQA1"}; $alleles{ "HLA_DQA1_".$fields[17] } = $genepos{"HLA_DQA1"};
    $alleles{ "HLA_DQB1_".$fields[18] } = $genepos{"HLA_DQB1"}; $alleles{ "HLA_DQB1_".$fields[19] } = $genepos{"HLA_DQB1"};
    $alleles{ "HLA_DRB1_".$fields[20] } = $genepos{"HLA_DRB1"}; $alleles{ "HLA_DRB1_".$fields[21] } = $genepos{"HLA_DRB1"};

    $alleles{ "HLA_A_".substr($fields[6],0,2) } = $genepos{"HLA_A"};     $alleles{ "HLA_A_".substr($fields[7],0,2) } = $genepos{"HLA_A"};
    $alleles{ "HLA_B_".substr($fields[8],0,2) } = $genepos{"HLA_B"};     $alleles{ "HLA_B_".substr($fields[9],0,2) } = $genepos{"HLA_B"};
    $alleles{ "HLA_C_".substr($fields[10],0,2) } = $genepos{"HLA_C"};    $alleles{ "HLA_C_".substr($fields[11],0,2) } = $genepos{"HLA_C"};
    $alleles{ "HLA_DPA1_".substr($fields[12],0,2) } = $genepos{"HLA_DPA1"}; $alleles{ "HLA_DPA1_".substr($fields[13],0,2) } = $genepos{"HLA_DPA1"};
    $alleles{ "HLA_DPB1_".substr($fields[14],0,2) } = $genepos{"HLA_DPB1"}; $alleles{ "HLA_DPB1_".substr($fields[15],0,2) } = $genepos{"HLA_DPB1"};
    $alleles{ "HLA_DQA1_".substr($fields[16],0,2) } = $genepos{"HLA_DQA1"}; $alleles{ "HLA_DQA1_".substr($fields[17],0,2) } = $genepos{"HLA_DQA1"};
    $alleles{ "HLA_DQB1_".substr($fields[18],0,2) } = $genepos{"HLA_DQB1"}; $alleles{ "HLA_DQB1_".substr($fields[19],0,2) } = $genepos{"HLA_DQB1"};
    $alleles{ "HLA_DRB1_".substr($fields[20],0,2) } = $genepos{"HLA_DRB1"}; $alleles{ "HLA_DRB1_".substr($fields[21],0,2) } = $genepos{"HLA_DRB1"};

    $nsamples++;
}
close(FILE);
#print STDERR "$nsamples individuals read\n";


my @all_alleles = sort {
    if (int($alleles{$a}) == int($alleles{$b})){
	return $a cmp $b;
    }else{
	return int($alleles{$a}) <=> int($alleles{$b});
    }
} keys %alleles;
my $al = "";

#Writes unique alleles to map file
open(F2,">$mapname") || die "can't open $mapname\n";
foreach my $allele ( @all_alleles ){
    @tmp = split('\_',$allele);
    $al=$tmp[2];
    if ($al ne "NA" && $al ne "" && $al ne "0" && $al ne "0 0") {
	print F2 "6\t" . $allele . "\t0\t" . $genepos{$tmp[0]."_".$tmp[1]} . "\n";
    }
}
close(F2);

open(FILE,"$pedname") || die "can't open $pedname\n";
while (<FILE>) {
    chomp;
    s/^\s+//;
    @fields = split(/\s+/);
    print "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]";

    &PrintGenotypes("A",$fields[6],$fields[7]);
    &PrintGenotypes("C",$fields[10],$fields[11]);
    &PrintGenotypes("B",$fields[8],$fields[9]);
    &PrintGenotypes("DRB1",$fields[20],$fields[21]);
    &PrintGenotypes("DQA1",$fields[16],$fields[17]);
    &PrintGenotypes("DQB1",$fields[18],$fields[19]);
    &PrintGenotypes("DPA1",$fields[12],$fields[13]);
    &PrintGenotypes("DPB1",$fields[14],$fields[15]);

    print "\n";
}
close(FILE);

sub hashValueAscendingNum {
    my @all_alleles = sort {$alleles{$a}.$a cmp $alleles{$b}.$b} keys %alleles;
}

use strict;
sub PrintGenotypes {
    my($locus, $allele1, $allele2) = @_;
    my $al1 = $allele1;
    my $al2 = $allele2;
    $allele1="HLA_".$locus."_".$allele1;
    $allele2="HLA_".$locus."_".$allele2;
    my $G1 = "0";
    my $G2 = "0";
    my @tmp = ();
    my $al = "";
    my $tmplocus = "";

    foreach my $allele ( @all_alleles ){
	my@tmp = split('\_',$allele);
	$al=$tmp[2];
	$tmplocus=$tmp[0]."_".$tmp[1];

        #if allele being checked is valid and is within specific locus
	if ($al ne "NA" && $al ne "" && $al ne "0" && $al ne "0 0" && $tmplocus eq "HLA_".$locus){
	    #if there is no missing types at a locus
	    if ($al1 ne "NA" && $al1 ne "" && $al1 ne "0" && $al2 ne "NA" && $al2 ne "" && $al2 ne "0"){
		#assign genotype for allele 1
		if ($al1 eq "NA" || $al1 eq "" || $al1 eq "0"){
		    $G1 = "0";
		}else{
		    #if alleles matches at 4 digit, or if first 2 digit of data matches allele
		    if ($al1 eq $al || substr($al1,0,2) eq $al){
			$G1 = "P";
		    #if data is 2 digit and allele is 4 and they match at 2, set as "unknown"
		    }elsif (length($al1) == 2 && length($al) == 4 && substr($al,0,2) eq $al1){
			$G1 = "0";
		    }else{
			$G1 = "A";
		    }
		}

		#assign genotype for allele 2
		if ($al2 eq "NA" || $al2 eq "" || $al2 eq "0"){
		    $G2 = "0";
		}else{
		    #if alleles matches at 4 digit, or if first 2 digit of data matches allele
		    if ($al2 eq $al || substr($al2,0,2) eq $al){
			$G2 = "P";
		    #if data is 2 digit and allele is 4 and they match at 2, set as "unknown"
		    }elsif (length($al2) == 2 && length($al) == 4 && substr($al,0,2) eq $al2){
			$G2 = "0";
		    }else{
			$G2 = "A";
		    }
		}
		if ($G1 eq "0" || $G2 eq "0"){
		    print "\t0 0";
		}else{
		    print "\t" . $G1 . " " . $G2;
		}
	    }else{
		#print 0's if alleles are not typed at locus
		print "\t0 0";
	    }
	}
    }
}
