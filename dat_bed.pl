#!/usr/bin/perl -w
use strict;
my $chr="";
open IN, "reference.fa.2.7.7.80.10.50.500.dat";
open OUT,">str_bed_file";
while(<IN>){
	if(/Sequence: (\d+)/){$chr=$1};
	if(/^\d+/){
		my @line=split(/\s+/,$_);
		##1.output chromosome, 2.start coordinete, 3.end coordinate, 4.period, 5.reference copy number, 9.str score, 15.repeat unit
		print OUT "$chr\t$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[11]\t$line[12]\t$line[13]\t$line[14]\n";
	}
}
close OUT;
close IN;
