#!/usr/bin/perl

###########################################################################################
#
# Created by Sherman Jia, 2012
# encodeVariants.pl
#
# This script encodes multi-allelic PLINK markers (amino acids and SNPs) into bi-allelic markers
# Input files include a normal PLINK .ped and .map file (where the .ped file contains
# multi-allelic positions).
#
###########################################################################################

if (scalar(@ARGV) !=3){ die "usage: ./encodeVariants.pl input.ped input.map output (.ped/.map)\n";}

$pedfile=shift;
$mapfile=shift;
$outputfile=shift;

my $nsamples = 0;
my $nsnps = 0;
my @alleles = ();
my @ID = ();
my @map = ();
my @SNP1 = ();
my @SNP2 = ();
my @multiallelic = ();
my %genepos = ();

open(FILE,"$mapfile") || die "can't open $mapfile\n";
while (my $c = <FILE>) {
    #Reads map file to get SNP names
    chomp($c);
    $nsnps++;
    $map[$nsnps] = $c;
}
close(FILE);
#print STDERR "$nsnps snps read\n";


open(FILE,"$pedfile") || die "can't open $pedfile\n";
while (<FILE>) {
    #Reads ped file to look for multi-allelic SNPs
    chomp;
    s/^\s+//;
    @fields = split(/\t/); 
    $nsamples++;
    $ID[$nsamples] = $fields[0];
    for (my $i = 1; $i < 6; $i++){
	$ID[$nsamples] = $ID[$nsamples] . "\t" . $fields[$i];
    }

    for (my $i = 6; $i < @fields; $i++){
	@SNP = split(/ /,$fields[$i]);
	if ($SNP[0] ne "0"){
	    if (rindex($alleles[$i-5],$SNP[0]) < 0 ){
		$alleles[$i-5] = $alleles[$i-5] . $SNP[0];
	    }
	}
	if ($SNP[1] ne "0"){
	    if (rindex($alleles[$i-5],$SNP[1]) < 0){
		$alleles[$i-5] = $alleles[$i-5] . $SNP[1];
	    }
	}
    }

}
close(FILE);


#Print new PED file
open(PEDOUT,">$outputfile".".ped") || die "can't open $outputfile.ped\n";
open(FILE,"$pedfile") || die "can't open $pedfile\n";
my $s = 0;
while (my $l=<FILE>) {
    #Reads ped file to look for multi-allelic SNPs
    chomp($l);
    @fields = split(/\t/,$l); 
    $s++;
    print PEDOUT $ID[$s];
    for (my $i = 6; $i < @fields; $i++){
	@SNP = split(/ /,$fields[$i]);
	if (length($alleles[$i-5]) > 2){
	    #Code each allele separately
	    for (my $j = 0; $j < length($alleles[$i-5]); $j++){
		if ($SNP[0] eq "0" || $SNP[1] eq "0"){
		    print PEDOUT "\t0 0";
		}else{
		    if (substr($alleles[$i-5],$j,1) eq $SNP[0]){
			print PEDOUT "\tP";
		    }else{
			print PEDOUT "\tA";
		    }
		    if (substr($alleles[$i-5],$j,1) eq $SNP[1]){
			print PEDOUT " P";
		    }else{
			print PEDOUT " A";
		    }
		}
	    }

	    if (length($alleles[$i-5]) > 3){
		#if > 3 alleles at site, create marker that tags pairs of alleles

		my $j_end = 0;
		if (length($alleles[$i-5]) == 4){
		    $j_end = 1;
		}else{
		    $j_end = length($alleles[$i-5]);
		}

		for (my $j = 0; $j < $j_end; $j++){
		    for (my $k = $j+1; $k < length($alleles[$i-5]); $k++){
			if ($SNP[0] eq "0" || $SNP[1] eq "0"){
			    print PEDOUT "\t0 0";
			}else{
			    if (substr($alleles[$i-5],$j,1) eq $SNP[0] || substr($alleles[$i-5],$k,1) eq $SNP[0]){
				print PEDOUT "\tP";
			    }else{
				print PEDOUT "\tA";
			    }
			    if (substr($alleles[$i-5],$j,1) eq $SNP[1] || substr($alleles[$i-5],$k,1) eq $SNP[1]){
				print PEDOUT " P";
			    }else{
				print PEDOUT " A";
			    }
			}
		    }
		}

		if (length($alleles[$i-5]) > 5){
		    #if > 5 alleles at site, create marker that tags 3 alleles

		    if (length($alleles[$i-5]) == 6){
			$j_end = 1;
		    }else{
			$j_end = length($alleles[$i-5]);
		    }

		    for (my $j = 0; $j < $j_end; $j++){
			for (my $k = $j+1; $k < length($alleles[$i-5]); $k++){
			    for (my $l = $k+1; $l < length($alleles[$i-5]); $l++){
				if ($SNP[0] eq "0" || $SNP[1] eq "0"){
				    print PEDOUT "\t0 0";
				}else{
				    if (substr($alleles[$i-5],$j,1) eq $SNP[0] || substr($alleles[$i-5],$k,1) eq $SNP[0] || substr($alleles[$i-5],$l,1) eq $SNP[0]){
					print PEDOUT "\tP";
				    }else{
					print PEDOUT "\tA";
				    }
				    if (substr($alleles[$i-5],$j,1) eq $SNP[1] || substr($alleles[$i-5],$k,1) eq $SNP[1] || substr($alleles[$i-5],$l,1) eq $SNP[1]){
					print PEDOUT " P";
				    }else{
					print PEDOUT " A";
				    }
				}
			    }
			}
		    }
		}
	    }
	}else{
	    if ($SNP[0] eq "0" || $SNP[1] eq "0"){
		print PEDOUT "\t0 0";
	    }else{
		print PEDOUT "\t" . $SNP[0] . " " . $SNP[1] ;
	    }
	}
    }
    print PEDOUT "\n";
}
close(FILE);
close(PEDOUT);

#Write new map file
open(MAPOUT,">$outputfile".".map") || die "can't open $outputfile.map\n";
for (my $i = 6; $i < @fields; $i++){
    if (length($alleles[$i-5]) > 2){
	my @fields = split(/\t/,$map[$i-5]);
	for (my $j = 0; $j < length($alleles[$i-5]); $j++){
	    print MAPOUT $fields[0] . "\t" . $fields[1] . "_" . substr($alleles[$i-5],$j,1) . "\t" . $fields[2] . "\t" . $fields[3] . "\n";
	}
	if (length($alleles[$i-5]) > 3){

	    my $j_end = 0;
	    if (length($alleles[$i-5]) == 4){
		$j_end = 1;
	    }else{
		$j_end = length($alleles[$i-5]);
	    }

	    for (my $j = 0; $j < $j_end; $j++){
		for (my $k = $j+1; $k < length($alleles[$i-5]); $k++){
		    print MAPOUT $fields[0] . "\t" . $fields[1] . "_" . substr($alleles[$i-5],$j,1).substr($alleles[$i-5],$k,1) . "\t" . $fields[2] . "\t" . $fields[3] . "\n";
		}
	    }
	    if (length($alleles[$i-5]) > 5){

		if (length($alleles[$i-5]) == 6){
		    $j_end = 1;
		}else{
		    $j_end = length($alleles[$i-5]);
		}

		for (my $j = 0; $j < $j_end; $j++){
		    for (my $k = $j+1; $k < length($alleles[$i-5]); $k++){
			for (my $l = $k+1; $l < length($alleles[$i-5]); $l++){
			    print MAPOUT $fields[0] . "\t" . $fields[1] . "_" . substr($alleles[$i-5],$j,1).substr($alleles[$i-5],$k,1).substr($alleles[$i-5],$l,1) . "\t" . $fields[2] . "\t" . $fields[3] . "\n";
			}
		    }
		}
	    }
	}
    }else{
	print MAPOUT $map[$i-5] . "\n";
    }
}
close(MAPOUT);
