#!/usr/bin/perl -w
use List::MoreUtils qw(uniq);
use List::Util qw/shuffle/;
use strict;

my $in="SSR.vcf";
my $ou="random";
my $n=3;
my $l=5;
my $m=5;

&parse_command_line();

my $line="";
my @col=();
my @ele=();
my @cor=();
my ($i,$j,$x,$y,$z)=("","","","","","");
my @tmp=();
my @u_tmp=();
my $TMP="";
my $per=0;

my %ind=();
my $rep="";
my @motif=();
my $c=0;
my @list=();
my @newlist=();
my @vol=();
my @u_vol=();
my @s_vol=();

open IN,$in or die "$!";
while($line=<IN>){
	if($line=~/^##/){next}
	elsif($line=~/^#CHROM/){
		@col=split(/\t/,$line);
		$per=scalar(@col)-9;
	}
	else{
		$c=$per;
		@tmp=();
		chomp $line;
		if($line=~/MOTIF=(\w+);/){
			$rep=$1;		
		}
		@col=split(/\t/,$line);
		@motif=split(/\,/,$col[4]);
		unshift(@motif,$col[3]);
		foreach $i (9..$#col){
			if($col[$i]=~/^\.\/\./ or length($col[$i])==0){
				$c--
			}
			else{
				@ele=split(/:/,$col[$i]);
				@cor=split(/\//,$ele[0]);
				push(@tmp,@cor);
			}
		}
		@u_tmp=uniq @tmp;
		if(length($rep)>=$n && $c==$per && scalar(@u_tmp)>=$l){
			print STDERR $col[0],":",$col[1],"\t@u_tmp\n";
			foreach $i (9..$#col){
				@ele=split(/:/,$col[$i]);
				@cor=split(/\//,$ele[0]);
				push(@{$ind{$col[0].":".$col[1]}{$i-8}},length($motif[$cor[0]]));
				push(@{$ind{$col[0].":".$col[1]}{$i-8}},length($motif[$cor[1]]));
			}
		}
	}
}
close IN;

$TMP=$ou."\.tmp";
open OUT,">$TMP";
	print OUT "Individual Number\tAllele Number\tRange\tMaximum\tMinimum\n";
open ALLELE,">Allele_number.txt";
	foreach $i ($m..$per-1){
		print ALLELE "S$i\t";
		}
	print ALLELE "S$per\n";
open RANGE,">Range.txt";
	foreach $i ($m..$per-1){
		print RANGE "S$i\t";
	}
	print RANGE "S$per\n";
open MAX,">Maximum.txt";
	foreach $i ($m..$per-1){
		print MAX "S$i\t";
	}
	print MAX "S$per\n";
open MINI,">Minimum.txt";
	foreach $i ($m..$per-1){
		print MINI "S$i\t";
	}
	print MINI "S$per\n";

for $i (1..$per){
	push(@list,$i);
}
foreach $i (keys %ind){
	print STDERR $i,"\n";
	my @newlist=shuffle @list;
	@vol=();
	foreach $j ($m..$per-1){
		foreach $z (0..$j-1){
			push(@vol,@{$ind{$i}{$newlist[$z]}});
		}
		@u_vol=uniq @vol;
		@s_vol=sort {$a<=>$b} @u_vol;
		print OUT $j,"\t",scalar(@s_vol),"\t",$s_vol[-1]-$s_vol[0],"\t",$s_vol[-1],"\t",$s_vol[0],"\n";
		print ALLELE scalar(@s_vol),"\t";
		print RANGE $s_vol[-1]-$s_vol[0],"\t";
		print MAX $s_vol[-1],"\t";
		print MINI $s_vol[0],"\t";
	}
	$j=$per;
	foreach $z (0..$j-1){
		push(@vol,@{$ind{$i}{$newlist[$z]}});
	}
	@u_vol=uniq @vol;
	@s_vol=sort {$a<=>$b} @u_vol;
	print OUT $j,"\t",scalar(@s_vol),"\t",$s_vol[-1]-$s_vol[0],"\t",$s_vol[-1],"\t",$s_vol[0],"\n";
	print ALLELE scalar(@s_vol),"\n";
	print RANGE $s_vol[-1]-$s_vol[0],"\n";
	print MAX $s_vol[-1],"\n";
	print MINI $s_vol[0],"\n";
}
close OUT;
close MINI;
close MAX;
close RANGE;
close ALLELE;

$TMP=$ou."\.plot\.R";
open OUT,">$TMP";
print OUT <<EOF;

#!/usr/bin/Rscript
library(ggplot2)
library("psych")
library("stats")
library("PMCMR")
tiff("$ou.pairs.panels.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
dt<-read.table("$ou.tmp",header=TRUE,sep="\t")
dt\$Individual.Number<-as.factor(dt\$Individual.Number)

pairs.panels(dt[c("Individual.Number","Allele.Number","Range","Maximum","Minimum")], stars=TRUE)
dev.off()
tiff("$ou.AlleleNumber_IndividualNumber.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
p<-ggplot(data=dt, aes(x=Individual.Number,y=Allele.Number,fill=Individual.Number))+geom_violin()+geom_jitter(alpha=0.3,col="red",show.legend=FALSE) +theme(legend.position='none')+xlab("Individual Number")+ylab("Allele Number")+theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"))
p
dev.off()
tiff("$ou.Range_IndividualNumber.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
p<-ggplot(data=dt, aes(x=Individual.Number,y=Range,fill=Individual.Number))+geom_violin()+geom_jitter(alpha=0.3,col="red",show.legend=FALSE) +theme(legend.position='none')+xlab("Individual Number")+ylab("Range")+theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"))
p
dev.off()
tiff("$ou.Maximum_IndividualNumber.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
p<-ggplot(data=dt, aes(x=Individual.Number,y=Maximum,fill=Individual.Number))+geom_violin()+geom_jitter(alpha=0.3,col="red",show.legend=FALSE) +theme(legend.position='none')+xlab("Individual Number")+ylab("Maximum")+theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"))
p
dev.off()
tiff("$ou.Minimum_IndividualNumber.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
p<-ggplot(data=dt, aes(x=Individual.Number,y=Minimum,fill=Individual.Number))+geom_violin()+geom_jitter(alpha=0.3,col="red",show.legend=FALSE) +theme(legend.position='none')+xlab("Individual Number")+ylab("Minimum")+theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"))
p
dev.off()

sink("$ou.test.txt")
Allele_number<-read.table("Allele_number.txt",header=TRUE,sep="\t")
Allele_number<-as.matrix(Allele_number)
friedman.test(Allele_number)
posthoc.friedman.nemenyi.test(Allele_number)
Range<-read.table("Range.txt",header=TRUE,sep="\t")
Range<-as.matrix(Range)
friedman.test(Range)
posthoc.friedman.nemenyi.test(Range)
Maximum<-read.table("Maximum.txt",header=TRUE,sep="\t")
Maximum<-as.matrix(Maximum)
friedman.test(Maximum)
posthoc.friedman.nemenyi.test(Maximum)
Minimum<-read.table("Minimum.txt",header=TRUE,sep="\t")
Minimum<-as.matrix(Minimum)
friedman.test(Minimum)
posthoc.friedman.nemenyi.test(Minimum)
sink()
EOF
close OUT;
system("chmod 755 $TMP");
system("R <$TMP --no-save >$ou.test.log");

sub parse_command_line
{
	if(scalar @ARGV ==0){
		&usage();
		exit(0)
	}
	else{
		while (@ARGV) {
			$_ = shift @ARGV;
			if( $_ =~ /^-i$/   ){ $in = shift @ARGV }			
			elsif( $_=~ /^-o$/ ){ $ou = shift @ARGV }
			elsif( $_=~ /^-n$/ ){ $n  = shift @ARGV }
			elsif( $_=~ /^-a$/ ){ $l  = shift @ARGV }
			elsif( $_=~ /^-m$/ ){ $m  = shift @ARGV }
			elsif( $_=~ /^-h$/ ){ &usage()		}
			else {
				print STDERR "Unknown command line option: '$_'\n";
				&usage();
				exit(0);
			}
		}
	}
}


sub usage
{
	print STDERR <<EOQ;

	random.pl
	
	version 1.1
	-i: the input file [SSR.vcf]
	-o: the prefix of output file [random]
	-n: the shortest length of motif [3]
	-a: miniumum number of allela [5]
	-m: minumum number of individual [5]

EOQ
}
