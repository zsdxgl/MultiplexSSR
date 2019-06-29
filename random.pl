#!/usr/bin/perl -w
use List::MoreUtils qw(uniq);
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
		if(length($rep)>$n && $c==$per && scalar(@u_tmp)>$l){
			print STDERR "@u_tmp\n";
			foreach $i (9..$#col){
				@ele=split(/:/,$col[$i]);
				@cor=split(/\//,$ele[0]);
				push(@{$ind{$col[0].":".$col[1]}{$i-9}},length($motif[$cor[0]]));
				push(@{$ind{$col[0].":".$col[1]}{$i-9}},length($motif[$cor[1]]));
			}
		}
	}
}
close IN;

$TMP=$ou."\.tmp";
open OUT,">$TMP";
print OUT "Individual Number\tAllele Number\tRange\tMaximum\tMinimum\n";
foreach  $i ($m..$per){
	foreach $j (keys %ind){
		print STDERR $i,"\t",$j,"\t";
		@vol=();
		@tmp=&Ra($i,$per);
		print STDERR "@tmp","\n";
		foreach $x (@tmp){
			push(@vol,@{$ind{$j}{$x}});
		}
		@u_vol=uniq @vol;
		@s_vol=sort {$a<=>$b} @u_vol;
		print OUT $i,"\t",scalar(@s_vol),"\t",$s_vol[-1]-$s_vol[0],"\t",$s_vol[-1],"\t",$s_vol[0],"\n";
	}
	print STDERR "\n";
}
close OUT;

$TMP=$ou."\.R";
open OUT,">$TMP";
print OUT <<EOF;

#!/usr/bin/Rscript
library("psych")
attach(airquality)
library(ggplot2)
tiff("$ou.pairs.panels.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
dt<-read.table("$ou.tmp",header=TRUE,sep="\t")
dt\$Individual.Number<-as.factor(dt\$Individual.Number)

pairs.panels(dt[c("Individual.Number","Allele.Number","Range","Maximum","Minimum")], stars=TRUE)
dev.off()
tiff("$ou.AlleleNumber.IndividualNumber.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
p<-ggplot(data=dt, aes(x=Individual.Number,y=Allele.Number,fill=Individual.Number))+geom_violin()+geom_jitter(alpha=0.3,col="red",show.legend=FALSE) +theme(legend.position='none')+xlab("Individual Number")+ylab("Allele Number")
p
dev.off()
tiff("$ou.Range.IndividualNumber.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
p<-ggplot(data=dt, aes(x=Individual.Number,y=Range,fill=Individual.Number))+geom_violin()+geom_jitter(alpha=0.3,col="red",show.legend=FALSE) +theme(legend.position='none')+xlab("Individual Number")+ylab("Range")
p
dev.off()
tiff("$ou.Maximum.IndividualNumber.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
p<-ggplot(data=dt, aes(x=Individual.Number,y=Maximum,fill=Individual.Number))+geom_violin()+geom_jitter(alpha=0.3,col="red",show.legend=FALSE) +theme(legend.position='none')+xlab("Individual Number")+ylab("Maximum")
p
dev.off()
tiff("$ou.Minimum.IndividualNumber.tiff",res=300,width = 4800, height = 4800, compression = "lzw")
p<-ggplot(data=dt, aes(x=Individual.Number,y=Minimum,fill=Individual.Number))+geom_violin()+geom_jitter(alpha=0.3,col="red",show.legend=FALSE) +theme(legend.position='none')+xlab("Individual Number")+ylab("Minimum")
p
dev.off()

sink("$ou.text.txt")
kruskal.test(Allele.Number~Individual.Number, data=dt)
kruskal.test(Range~Individual.Number, data=dt)
kruskal.test(Maximum~Individual.Number, data=dt)
kruskal.test(Minimum~Individual.Number, data=dt)
pairwise.wilcox.test(dt\$Allele.Number,dt\$Individual.Number,p.adj = "bonf")
pairwise.wilcox.test(dt\$Range,dt\$Individual.Number,p.adj = "bonf")
pairwise.wilcox.test(dt\$Maximum,dt\$Individual.Number,p.adj = "bonf")
pairwise.wilcox.test(dt\$Minimum,dt\$Individual.Number,p.adj = "bonf")
sink()
EOF
close OUT;
system("chmod 755 $ou\.R");
system("R <$ou\.R --no-save >$ou.test.log");


sub Ra{
	my ($i,$j)=@_;
	my %sns = ();
	my $TMP="";
	for(my $z = 0;$z < 100000; $z++) {
		$TMP = int(rand($j));
		$sns{$TMP}=1;
		my @x=keys %sns;
		if(scalar(@x)==$i){
			return(@x);
			last
		}
	}
}

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
	-n: the short length of motif [3]
	-a: miniumum number of allela [5]
	-m: minumum number of individual [5]

EOQ
}
