#!/usr/bin/perl -w
use strict;
use Benchmark;
use List::Util qw(max min);


########path
my $bin="";
my $trf="";
my $lobSTR="";
my $primer3="";
my $ePCR="";
my $Multiplx="";
########
my $bwa=`which bwa`;
chomp $bwa;
my $samtools=`which samtools`;
chomp $samtools;
my $bcftools=`which bcftools`;
chomp $bcftools;

###reference assembly
my $ref="./reference.fa";
###Samples
my $samples="./rawdata";
my $list="./input.txt";
####
my $t=4;
my $T=60.0;
my $p="Multiplex";
my $l=3;
my $n=5;
my $al=5;
my $d=1;
my $mi=80;
my $ma=470;
my $g=10;
my $s="F";

&parse_command_line();

####
$trf=$bin."/trf409.linux64";
$lobSTR=$bin."/lobSTR";
$primer3=$bin."/primer3-2.4.0/src/primer3_core";
$ePCR=$bin."/ePCR";
$Multiplx=$bin."/Multiplx";
####
my $i="";
my $line="";
my @file=();
my $tmp="";
my $tmpfile="";
my $str="";
my @TMPFILE=();
my @tmpelement=();
my $cmd="";
my %ssr=();

##################################
if($s eq "T"){
	goto LABEL
}
###
#goto SSR;
$tmp='./SNP';
if(-e $tmp){
	$cmd='rm -rf '.$tmp;
	system($cmd);
	$cmd='mkdir SNP';
	system($cmd);
}else{
	$cmd='mkdir SNP';
	system($cmd);
};
$tmp='./SSR';
if(-e $tmp){
	$cmd='rm -rf '.$tmp;
	system($cmd);
	$cmd='mkdir SSR';
	system($cmd);
}else{
	$cmd='mkdir SSR';
	system($cmd)
}


###########
###Calling SNPs
##########
unless((-e $ref.".amb") and (-e $ref.".ann") and (-e $ref."bwt") and (-e $ref."sa") ){	
$cmd=$bwa.' index '.$ref;
system($cmd);
}
open IN, $list or die "$!";
while($line=<IN>){
	chomp $line;
	@file=split(/\t/,$line);
	($tmp)=$file[0]=~/(\S+)[\._]1\.\w+/;
	$cmd=$bwa.' mem -t '.$t.' '.$ref.' '.$samples.'/'.$file[0].' '.$samples.'/'.$file[1].'|'.$samtools.' view -b -f 0x2 -F 0x100 -F 0x800 |'.$samtools.' sort --threads '.$t.' > ./SNP/'.$tmp.'.SNP.bam';
	system($cmd);
}
$cmd='ls '.'./SNP/*.SNP.bam >bam.list';
print STDERR $cmd,"\n";
system($cmd);
$cmd=$bcftools.' mpileup -C 50 -f '.$ref.' -b bam.list -q 20 -Q 20 -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR | '.$bcftools.' call -cv - > SNP.vcf';
print STDERR $cmd,"\n";
system($cmd);

#SSR:
####################################
#######mask reference and detect SSR
####################################
$cmd=$trf.' '.$ref.' 2 7 7 80 10 50 500 -f -d -m';
system($cmd);
$cmd='rm *html';
system($cmd);

&dat_bed();

##################################
####################SSR genotyping
##################################
##Build the lobSTR index
$cmd='python '.$lobSTR.'/share/lobSTR/scripts/lobstr_index.py --str ./SSR/str_bed_file --ref '.$ref.' --out ./SSR/index';
print STDERR $cmd,"\n";
system($cmd);
#Generating the STRInfo file
$cmd='python '.$lobSTR.'/share/lobSTR/scripts/GetSTRInfo.py ./SSR/str_bed_file '.$ref.' > ./SSR/strinfo.tab';
print STDERR $cmd,"\n";
system($cmd);
#genotyping SSR
open IN, $list or die "$!";
while($line=<IN>){
	chomp $line;
	@file=split(/\t/,$line);
	($tmp)=$file[0]=~/(.*)[\._]1\.\w+/;
	$cmd=$lobSTR.'/bin/lobSTR --index-prefix ./SSR/index/lobSTR_   --p1 '.$samples.'/'.$file[0].' --p2 '.$samples.'/'.$file[1].' -q -p '.$t.' -u  --rg-sample '.$tmp.' --rg-lib '.$tmp.' --out ./SSR/'.$tmp.'.SSR';
	print STDERR $cmd,"\n";
	system($cmd);
	$cmd=$samtools.' sort ./SSR/'.$tmp.'.SSR.aligned.bam -T tmp -o ./SSR/'.$tmp.'.SSR.sorted.bam';
	print STDERR $cmd,"\n";
	system($cmd);
	$cmd=$samtools.' index ./SSR/'.$tmp.'.SSR.sorted.bam';
	print STDERR $cmd,"\n";
	system($cmd);
	$tmpfile='./SSR/'.$tmp.'.SSR.sorted.bam';
	push(@TMPFILE,$tmpfile);
}
$str=join(',',@TMPFILE);
$cmd=$lobSTR.'/bin/allelotype --command classify --bam  '.$str.' --index-prefix ./SSR/index/lobSTR_ --strinfo ./SSR/strinfo.tab --noise_model '.$lobSTR.'/share/lobSTR/models/illumina_v3.pcrfree --min-het-freq 0.2  --min-border 5 --min-bp-before-indel 7 --maximal-end-match 15 --min-read-end-match 5 --max-matedist 1000 --filter-clipped --noweb --out SSR';
print STDERR $cmd,"\n";
system($cmd);
################################
LABEL:
################################
###############Desinging primers
################################
print STDERR "extract_flank()\n";
&extract_flank();
print STDERR "DesignPrimer()\n";
&DesignPrimer();
print STDERR "formatchange()\n";
&formatchange();
###Filtering primers
$cmd=$ePCR.'/famap -tn -b genome.famap '.$ref;
print STDERR $cmd,"\n";
system($cmd);
$cmd=$ePCR.'/fahash -b genome.hash -w 12 -f 3 ./genome.famap';
print STDERR $cmd,"\n";
system($cmd);
$cmd=$ePCR.'/re-PCR -S genome.hash -o ePCR.out.txt -n 3 -g 3 -d '.$mi.'-'.$ma.' -r + primer.input4ePCR.txt';
print STDERR $cmd,"\n";
system($cmd);
&extract();

################################
################grouping primers
################################
$cmd='ln -s '.$Multiplx.'/thermodynamics.txt ./';
system($cmd);
$cmd=$Multiplx.'/cmultiplx -primers primer.ff.txt -calcscores 12345 -saveprimscores primer-scores.txt -saveprodscores product-scores.txt';
print STDERR $cmd,"\n";
system($cmd);
$cmd=$Multiplx.'/cmultiplx -primers primer.ff.txt -loadprimscores primer-scores.txt -loadprodscores product-scores.txt -initialorder friends -stringency normal  -calcgroups 50 1000 -savegroups groups.txt';
print STDERR "\n";
system($cmd);
$cmd='dos2unix groups.txt';
print STDERR $cmd,"\n";
system($cmd);
&select();
##################################
$cmd='rm thermodynamics.txt';
system($cmd);
#$cmd='rm ePCR.out.txt';
#system($cmd);
#$cmd='rm bam.list';
#system($cmd);
#$cmd='rm flanking_seq.txt';
#system($cmd);
#$cmd='rm genome.famap';
#system($cmd);
#$cmd='rm genome.hash';
#system($cmd);
#$cmd='rm groups.txt';
#system($cmd);
#$cmd='rm primer.f.txt';
#system($cmd);
#$cmd='rm primer.input4ePCR.txt';
#system($cmd);
#$cmd='rm primer-scores.txt';
#system($cmd);
#$cmd='rm product-scores.txt';
#system($cmd);
#$cmd='rm reference.fa.2.7.7.80.10.50.500.dat';
#system($cmd);
#$cmd='rm reference.fa.2.7.7.80.10.50.500.mask';
#system($cmd);
#$cmd='rm reference.fa.amb';
#system($cmd);
#$cmd='rm reference.fa.ann';
#system($cmd);
#$cmd='rm reference.fa.bwt';
#system($cmd);
#$cmd='rm reference.fa.fai';
#system($cmd);
#$cmd='rm reference.fa.flat';
#system($cmd);
#$cmd='rm reference.fa.gdx';
#system($cmd);
#$cmd='rm reference.fa.pac';
#system($cmd);
#$cmd='rm reference.fa.sa';
#system($cmd);
#$cmd='rm SSR.allelotype.stats';
#system($cmd);
############################################
sub dat_bed
{
my $chr="";
	$tmp=$ref.'.2.7.7.80.10.50.500.dat';
	open IN,$tmp;
	open OUT,">./SSR/str_bed_file";
	while(<IN>){
		chomp;
	        if(/Sequence: (.*)$/){$chr=$1};
	        if(/^\d+/){
                	my @line=split(/\s+/,$_);
                	##1.output chromosome, 2.start coordinete, 3.end coordinate, 4.period, 5.reference copy number, 9.str score, 15.repeat unit
               		print OUT "$chr\t$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[11]\t$line[12]\t$line[13]\t$line[14]\n";
        	}
	}
	close OUT;
	close IN;
}
##############################################
sub extract_flank
{


	my $out  ="flanking_seq.txt";
	my $SNP  ="SNP.vcf";
	my $SSR  ="SSR.vcf";
	my $refM  =$ref.'.2.7.7.80.10.50.500.mask';


	########################
	my $line;my @vol;my @cell=();my @ele=();
	my $i="";my @allele=();my @uniqallele=();
	my $tmp="";my %seen=();
	my $ind=0;
	my %REF=();my $chr="";
	open REF, "$refM" or die "$!";
        	while($line=<REF>){
                	chomp $line;
                	if($line=~/^>(.+)$/){
                        	if(length($tmp)>0 and length($chr)>0){
                                	print STDERR length($tmp),"\n";
                                	@vol=split(//,$tmp);
                                	foreach $i (0..$#vol){
                                       		$REF{$chr}{$i+1}=$vol[$i];
                                        	# print STDERR $REF{$chr}{$i+1};
                                	}
                        	}
                        	$chr=$1;
                        	print STDERR $chr,"\t";
                        	$tmp="";
                	}
                	else{
                        	$tmp.=$line;
                	}
        	}
        	if(length($tmp)>0 and length($chr)>0){
                	print STDERR length($tmp),"\n";
                	@vol=split(//,$tmp);
                	foreach $i (0..$#vol){
                        	$REF{$chr}{$i+1}=$vol[$i];
                        	#print STDERR $REF{$chr}{$i+1};
                	}
        	}
	close REF;
	#################
	#input SNP.vcf and mask mutation
	open M,"$SNP" or die "$!";
        	while($line=<M>){
                	chomp $line;
                	if($line=~/^#/){next}
                	else{
                        	@vol=split(/\t/,$line);
                        	if($vol[5]>20){
                                	$REF{$vol[0]}{$vol[1]}="N";
                        	}
                	}
        	}
	close M;
	##################
	#extract flanking sequences
	open SSR,$SSR or die "$!";
	open OUT, ">$out";
			#print OUT "chr;start-end;MOTIF;REF;alleles\n";
        	while($line=<SSR>){
                	chomp $line;
					@allele=();
					%seen=();
					@cell=();
					@ele=();
					$ind=0;
                	if($line=~/^#/){next}
                	else{
                        	@vol=split(/\t/,$line);
                        	## QUAL, motif length, sample number with data, number of allele for selection
                        	## start and end of variance, length of ref, alleles for record
                        	#END=169898;MOTIF=TG;NS=10;REF=14;RL=28;RU=TG;VT=STR;RPA=10,12,13
                        	if($line=~/END=(\d+);MOTIF=(\w+);NS=(\d+);REF=(\d+);RL=(\d+);RU=\w+;VT=STR;RPA=(\d+.+)\tGT/){
                                	if(length($2)>=$l){
                                		for $i (9..$#vol){
							if($vol[$i]=~/:/){
								@cell=split(/:/,$vol[$i]);
								if($cell[4]>=$d){
									@ele=split(/\//,$cell[0]);
									push(@allele,@ele);
									$ind++;
								}
							}
						}
						foreach $i (@allele){
							$seen{$i}++
						}
						@uniqallele=sort keys %seen;
						print STDERR $vol[5],"\t",$2,"\t",$ind,"\t","@uniqallele","\t";
						if( defined($REF{$vol[0]}{$vol[1]-500})){print STDERR "-500","\t"}else{"N\t"}
						if(defined($REF{$vol[0]}{$1+500})){print STDERR "+500\t"}else{"N\t"};
                                		if($vol[5]>50 && length($2)>=$l && $ind >=$n && scalar(@uniqallele)>=$al  && defined($REF{$vol[0]}{$vol[1]-500}) && defined($REF{$vol[0]}{$1+500})){
                                       	 		print OUT ">$vol[0];$vol[1]-$1;$2;$4,@uniqallele\tStartpos_in_parent=501 Startpos_here=501 Length=$5\n";
                                       	 		foreach $i ($vol[1]-500 .. $1+500){
                                                		print OUT $REF{$vol[0]}{$i}
							}
                                        		print OUT "\n";
                                        		print STDERR "selected\n";
                                		}else{
                                        		print STDERR "\n";
                                		}
					}
                        	}
                	}
        	}
	close OUT;
	close SSR;
}

################################
sub formatchange
{
	my $line;my @vol=();my $ID="";
	$/="\n";
	open IN, "primer.f.txt" or die "$!";
	open OUT,">primer.input4ePCR.txt" or die "$!";
		while($line=<IN>){
		        chomp $line;
	       		@vol=split(/\t/,$line);
        		($ID)=$vol[0]=~/^(\d*)_L/;
			print STDERR "$ID\t$vol[1]\t$vol[3]\t$vol[9]\n";
        		print OUT "$ID\t$vol[1]\t$vol[3]\n";
		}
	close OUT;
	close IN;
}

sub extract
{
	my $line="";my %count=();my %sy=();my %start=();my %end=();my %pos=();
	my %len=();my @vol=();my $tmpL=0;my $i="";
	open EPCR,"<ePCR.out.txt";
	while($line=<EPCR>){
       		chomp $line;
        	if($line=~/^#/){next}
        	else{
               		@vol=split(/\t/,$line);
                	$count{$vol[0]}++;
                	$sy{$vol[0]}=$vol[2];
                	$pos{$vol[0]}=$vol[1];
                	$start{$vol[0]}=$vol[3];
                	$end{$vol[0]}=$vol[4];
                	($tmpL)=($line=~/\t(\d+)\/\d+\-\d+$/);
#               	print $tmpL,"\n";
                	$len{$vol[0]}=$tmpL;
        	}
	}
	close EPCR;

	my $ID;my %primer=();
	open PRIMER,"primer.f.txt";
	while($line=<PRIMER>){
        	chomp $line;
        	@vol=split(/\t/,$line);
        	($ID)=($line=~/^(\d+)\_/);
        	$primer{$ID}=$line;
	}
	close PRIMER;
	##ref
	my %REF=();my $chr="";
	open REF, $ref or die "$!";
		while($line=<REF>){
     			chomp $line;
     			if($line=~/>(.*)$/){$chr=$1}
     			else{
        			$REF{$chr}.=$line;
     			}
  		}
	close REF;
	my $M13="";
	open OUT,">primer.ff.txt";
	foreach  $i (sort {$a cmp $b} keys %count){
        	if($count{$i}==1 && $len{$i}>=$mi && $len{$i}<=$ma && $sy{$i} eq "+" ){
                	@vol=split(/\t/,$primer{$i});
                	#print OUT "$vol[0]\t$M13$vol[1]\t$vol[2]\t$M13$vol[3]";
                	#foreach my $j (4 .. $#vol){print OUT "\t$vol[$j]"}
                	#print OUT "\n";
                	print OUT $i,"\t",$vol[1],"\t",$vol[3],"\t",substr($REF{$pos{$i}},$start{$i}-1,$end{$i}-$start{$i}+1);
                	print OUT "\n";
        	}
	}
	close OUT;
}


sub DesignPrimer
{

	#Get Settings
	my ($path,$input,$output,$prefix,$appendForward,$leftTag,$rightTag,$detailLevel,$productSize,$optSize,$minSize,$maxSize,$optTM,$minTM,$maxTM,$diffTM,$optGC,$minGC,$maxGC,$maxPoly) = GetSettings();
	unless(-e $input)
	{
		print("File containing the input sequences could not be found. DesignPrimer aborted\n");
		print("Edit the file DesignPrimerInput.txt to specify an input file\n");
		print("Also make sure that this file is in the same directory as the script\n");
	}
	unless(-e $path)
	{
		print("Primer 3 could not be found. Process aborted\n");
                print("Edit the file DesignPrimerInput.txt to specify the position of Primer 3\n");
                print("Also make sure that this file is in the same directory as the script\n");
	}
	unless($productSize && $optSize && $minSize && $maxSize && $optTM && $minTM && $maxTM && $diffTM && $optGC && $minGC && $maxGC && $maxPoly)
	{
        print("Some parameters for Primer 3 are not valid\n");
        print "Edit the file DesignPrimerInput.txt and specify valid parameters for Primer3\n";
        print "if uncertain check the Primer 3 manual\n";
	}

#open fasta File
my($inputFH);
open($inputFH,$input);

#open outputFile
my $OUTPUT;
open($OUTPUT,">$output");

my $count=1;

	my(@cont);
        while(@cont = GetNextFastaRecord($inputFH))
        {
                my($fastaID,$dna)=@cont;
                #If the records are not defined cancel the procedure;
                if(!defined($fastaID) || !defined($dna)) {last;}

                #Retrieve the microsatellite target information from the fasta ID
                my($targetStart,$targetLength)=ParseSciRoKoFastaID($fastaID);

                #Start primer3 and retrieve the information
                @cont=GetPrimer3Results($dna,$targetStart,$targetLength,$path,$productSize,$optSize,$minSize,$maxSize,$optTM,$minTM,$maxTM,$diffTM,$optGC,$minGC,$maxGC,$maxPoly);
                #print @cont,"\n";
                my($fwdPrim,$revPrim,$fwdTM,$revTM,$prodSize)=ParsePrimer3Results(@cont);
                #print @cont,"\n";
                if(defined($fwdPrim) && defined($revPrim)){
                    print $OUTPUT $prefix.$count.$leftTag."\t".$appendForward.$fwdPrim."\t".$prefix.$count.$rightTag."\t".$revPrim;

                    if($detailLevel){
                        print $OUTPUT "\tTM_left:\t$fwdTM\tTM_right:\t$revTM\tPCR_product_size:\t$prodSize\tFull_fasta_ID:\t$fastaID";
                    }

                    print $OUTPUT "\n";
                        $count++;
                }
                else{
                    print "Could not find primer for: $fastaID\n";
                }
        }

	close($inputFH);
	close($OUTPUT);
	$cmd='rm temp.txt';
	system($cmd);
	$cmd='rm file';
	system($cmd);

	$count--;
	print "Successfully finished primer design\n";
	print "$count primer pairs have been designed\n";
}

sub ParsePrimer3Results
{
my(@content)=@_;
my($fwdPrim,$revPrim,$fwdTM,$revTM,$prodSize);
my $cont=join("\n",@content);

        ($fwdPrim)=$cont=~m/PRIMER_LEFT_\d*_SEQUENCE=(\S+)\s*\n/s;
        ($revPrim)=$cont=~m/PRIMER_RIGHT_\d*_SEQUENCE=(\S+)\s*\n/s;
        ($fwdTM)=$cont=~m/PRIMER_LEFT_\d*_TM=(\S+)\s*\n/s;
        ($revTM)=$cont=~m/PRIMER_RIGHT_\d*_TM=(\S+)\s*\n/s;
        ($prodSize)=$cont=~m/PRIMER_PAIR_\d*_PRODUCT_SIZE=(\S+)\s*\n/s;
    return ($fwdPrim,$revPrim,$fwdTM,$revTM,$prodSize);
}

sub GetPrimer3Results
{
        my($dna,$targetStart,$targetLength,$path,$productSize,$optSize,$minSize,$maxSize,$optTM,$minTM,$maxTM,$diffTM,$optGC,$minGC,$maxGC,$maxPoly)=@_;
    my $temp="temp.txt";
    my $OUTPUT;
    open($OUTPUT,">$temp");
    print $OUTPUT "SEQUENCE_ID=temp\n";
    print $OUTPUT "SEQUENCE_TEMPLATE=$dna\n";
    print $OUTPUT "SEQUENCE_TARGET=$targetStart,$targetLength\n";
    print $OUTPUT "PRIMER_PRODUCT_SIZE_RANGE=$productSize\n";
    print $OUTPUT "PRIMER_OPT_SIZE=$optSize\n";
    print $OUTPUT "PRIMER_MIN_SIZE=$minSize\n";
    print $OUTPUT "PRIMER_MAX_SIZE=$maxSize\n";
    print $OUTPUT "PRIMER_OPT_TM=$optTM\n";
    print $OUTPUT "PRIMER_MIN_TM=$minTM\n";
    print $OUTPUT "PRIMER_MAX_TM=$maxTM\n";
    print $OUTPUT "PRIMER_MAX_DIFF_TM=$diffTM\n";
    print $OUTPUT "PRIMER_OPT_GC_PERCENT=$optGC\n";
    print $OUTPUT "PRIMER_MIN_GC_PERCENT=$minGC\n";
    print $OUTPUT "PRIMER_MAX_GC_PERCENT=$maxGC\n";
    print $OUTPUT "PRIMER_MAX_POLY_X=$maxPoly\n";
    print $OUTPUT "PRIMER_NUM_RETURN=1\n";
    print $OUTPUT "=\n";
    close($OUTPUT);
    my $out="file";
    system (" $path  -default_version=1 < $temp >$out");
    my @content;
    open(IN, "$out") or die "$!";
    while (<IN>) {
        push @content, $_;
    }
    close IN;
    return @content;
    #print @content,"\n";

}

sub ParseSciRoKoFastaID
{
        my($fastaID)=@_;
        my($targetStart,$targetLength);
        ($targetStart)=$fastaID=~/Startpos_here=(\S*)/;  #Startpos_here=201  Length=20
        ($targetLength)=$fastaID=~/Length=(\S*)/;

        return($targetStart,$targetLength);
}

sub GetNextFastaRecord
{
 use strict;
        my($HANDLE)=@_;
        my ($defaultEnd)=$/;
        $/=">";
        my $line;

        do
        {
           $line=<$HANDLE>;
           unless($line){return;}

        }
        while($line=~/^\s*>\s*$/s);

        my($id,$sequence)=$line=~ /^(.*)\n((.*\n)*.*)/m;
        $sequence=~s/[\s>]//g;

        $/=$defaultEnd;
        return ($id,$sequence);

}

sub GetSettings
{
        my ($path,$input,$output,$prefix,$appendForward,$productSize);
        my ($optSize,$minSize,$maxSize,$optTM,$minTM,$maxTM,$diffTM,$optGC,$minGC,$maxGC);
        my ($maxPoly,$leftTag,$rightTag,$detailLevel);

        $path=$primer3;
        $input="flanking_seq.txt";
        $output="primer.f.txt";
        $prefix="";
        $appendForward="";
        $productSize="$mi-$ma";
       	$optSize=22;
        $minSize=18;
        $maxSize=26;
        $optTM=$T;
        $minTM=$T-3;
        $maxTM=$T+3;
        $optGC=50.0;
        $minGC=40.0;
        $maxGC=60.0;
        $maxPoly=3;
        $diffTM=3;
        $leftTag="_L";
        $rightTag="_R";
        $detailLevel=1;

    return ($path,$input,$output,$prefix,$appendForward,$leftTag,$rightTag,$detailLevel,$productSize,$optSize,$minSize,$maxSize,$optTM,$minTM,$maxTM,$diffTM,$optGC,$minGC,$maxGC,$maxPoly);
}

sub select
{
	my $line="";my @vol=();my %group=();
	%ssr=();my $id="";my $chr="";my $pos="";my $unit="";my $Rperiod="";my $Aperiod="";
	my @Ap=();my $max="";my $min="";

	open O,"primer.f.txt"  or die "$!";
        	while($line=<O>){
                	chomp $line;
                	@vol=split(/[\t\r]/,$line);
                	($id)=$line=~/^(\d+)_L/;
                	($chr,$pos,$unit,$Rperiod,$Aperiod)=$vol[11]=~/(\d+)\;(\d+-\d+)\;(\w+)\;(\d+),(.+)$/;
                	@Ap=split(/\s/,$Aperiod);
                	push (@Ap,$Rperiod);
                	$max= max @Ap;
                	$min= min @Ap;
                	$ssr{$id}{"unit"}=$unit;
                	@{$ssr{$id}{"period"}}=@Ap;
                	$ssr{$id}{"chr"}=$chr;
                	$ssr{$id}{"F_seq"}=$vol[1];
                	$ssr{$id}{"R_seq"}=$vol[3];
                	$ssr{$id}{"F_T"}=$vol[5];
                	$ssr{$id}{"R_T"}=$vol[7];
                	$ssr{$id}{"short"}=sprintf("%.0f",$vol[9]-(($Rperiod-$min)*length($unit)));
                	#print STDERR $id,"\t",$ssr{$id}{"short"} ,"\n";
                	$ssr{$id}{"long"}=sprintf("%.0f",$vol[9]+(($max-$Rperiod)*length($unit)));
        	}
	close O;

	open S,"primer.ff.txt" or die $!;
        	while($line=<S>){
                	chomp $line;
                	@vol=split(/\t/,$line);
                	$ssr{$vol[0]}{"seq"}=$vol[3];
        	}
	close S;

	open G,"groups.txt" or die "$!";
		$tmp=$p.'.summary.txt';
        	unlink($tmp);
		open SUM,">$tmp";
			print SUM "###FORMAT\tGroup name(Group Number)\tNumber of unique contig\n###the range of the SSR\n";
		close SUM;
		$tmp=$p.'.primer.txt';
		unlink($tmp);
		open DE,">$tmp";
			print DE "###SSR\trepeat unit\talleles\tforward primer\tTm\treverse primer\tTm\tthe seuqence of amplicon\n";
		close DE;
        	while($line=<G>){
                	chomp $line;
                	if($line=~/^Group/){
				@tmpelement=();
                        	@vol=split(/\t/,$line);
				foreach my $i (1..$#vol){
					push(@tmpelement,$vol[$i]);
				}
				if($#vol>2){
                        		&Group($vol[0],@tmpelement);
				}
                	}
        	}
	close G;
}

sub Group
{
        my $group="",my @GROUP; my @CHR=();my @tmpCHR=();my $chrN=0;my %count=();
        my $i="";my $j="";my $l=0;
	$tmp=$p.'.summary.txt';
        open SUM,">>$tmp";
	$tmp=$p.'.primer.txt';
        open DE,">>$tmp";
        ($group,@{$GROUP[0]})=@_;
	print STDERR $group,"\t","@{$GROUP[0]}","\n";
        @{$GROUP[0]} = sort {$ssr{$a}{"short"} <=> $ssr{$b}{"short"}} @{$GROUP[0]};
        foreach $i (@{$GROUP[0]}){push @CHR, $ssr{$i}{"chr"}};
        #print STDERR "@CHR","\n";
        @tmpCHR = grep { ++$count{ $_ } < 2; } @CHR;
        $chrN=@tmpCHR;
        print SUM $group,"\t",$chrN,"\n";
        print DE $group,"\t",$chrN,"\n";
        foreach $i (0 .. 30){
                if(defined ${$GROUP[$i]}[0]){
                        $l=@{$GROUP[$i]};
                        foreach $j (1 .. ($l-1)){
                                if(defined($GROUP[$i][$j]) && $ssr{$GROUP[$i][$j]}{"short"} < ($ssr{$GROUP[$i][$j-1]}{"long"} + $g)){
                                        #print $i,"\t",$j,"\t",$GROUP[$i][$j],"\t",$l,"\n";
                                        push @{$GROUP[$i+1]},$GROUP[$i][$j];
                                        splice (@{$GROUP[$i]},$j,1);
                                        $l=$l-1;
                                        redo;
    
                                 }
                        }
                }
                else{
                        last;
                }
        }
        foreach $i (0 .. 30){
                if(defined $GROUP[$i][0]){
                        print SUM "\t",'G',$i,"\t","@{$GROUP[$i]}","\n";
                        print SUM "\t\t";
                        foreach $j (0 .. $#{$GROUP[$i]}){
                                print SUM $ssr{$GROUP[$i][$j]}{"short"},"-",$ssr{$GROUP[$i][$j]}{"long"},"\t";
                                print DE $GROUP[$i][$j],"\t",$ssr{$GROUP[$i][$j]}{"chr"},"\t",$ssr{$GROUP[$i][$j]}{"unit"},"\t","@{$ssr{$GROUP[$i][$j]}{'period'}}","\t",$ssr{$GROUP[$i][$j]}{"F_seq"},"\t",$ssr{$GROUP[$i][$j]}{"F_T"},"\t",$ssr{$GROUP[$i][$j]}{"R_seq"},"\t",$ssr{$GROUP[$i][$j]}{"R_T"},"\t",$ssr{$GROUP[$i][$j]}{"seq"},"\n";
                        }
                        print SUM "\n";
                }
                else{
                        last;
                }
        }
        close SUM;
        close DE;
}


sub parse_command_line 
{
    if(scalar @ARGV ==0){
			usage();
			exit(0)
	}
    else{
        while (@ARGV) {
            $_ = shift @ARGV;
			if($_ =~ /^--bin$/)			{ $bin		= shift @ARGV; }
			elsif ($_ =~ /^--bwa$/) 	{ $bwa		= shift @ARGV; }
			elsif ($_ =~ /^--samtools$/){ $samtools	= shift @ARGV; }
			elsif ($_ =~ /^--bcftools$/){ $bcftools	= shift @ARGV; }
			elsif ($_ =~ /^--ref$/)		{ $ref		= shift @ARGV; }
			elsif ($_ =~ /^--samples$/)	{ $samples	= shift @ARGV; }
			elsif ($_ =~ /^--list$/)	{ $list		= shift @ARGV; }
			elsif ($_ =~ /^-t$/)		{ $t		= shift @ARGV; }
			elsif ($_ =~ /^-T$/)		{ $T		= shift @ARGV; }
			elsif ($_ =~ /^-p$/)		{ $p		= shift @ARGV; }
			elsif ($_ =~ /^-l$/)		{ $l		= shift @ARGV; }
			elsif ($_ =~ /^-n$/)		{ $n		= shift @ARGV; }
			elsif ($_ =~ /^-a$/)		{ $al		= shift @ARGV; }
			elsif ($_ =~ /^-d$/)		{ $d		= shift @ARGV; }
			elsif ($_ =~ /^-mi$/)		{ $mi		= shift @ARGV; }
			elsif ($_ =~ /^-ma$/)		{ $ma		= shift @ARGV; }
			elsif ($_ =~ /^-g$/)		{ $g		= shift @ARGV; }
			elsif ($_ =~ /^-s$/)		{ $s		= shift @ARGV; }
			elsif ($_ =~ /^-h$/)		{
								&usage()
			}
            else {
				print STDERR "Unknown command line option: '$_'\n";
				&usage();exit(0);
			}
		}
	}
	if(length($bin)==0 or length($samples)==0 or length($ref)==0 or length($list)==0){
		print STDERR "The path must be assigned for --bin, --samples, --ref, and --list\n";
		&usage();
		exit(0);
	}
}

sub usage 
{
    print STDERR <<EOQ;

	MultiplexSSR.pl

	version 1.1
	Path for softwares:
	--bin: path to the directory including lobSTR, trf, primer3_core, e-PCR, fafash, famap, re-PCR, and cmultiplx [\$PATH/MultiplexSSR/bin]
	--bwa: path to the program <bwa> [/usr/bin/bwa]
	--samtools: path to the program <samtools> [/usr/bin/samtools]
	--bctools: path to the program <bcftools> [/usr/bin/bcftools]

	Path for the reference
	--ref: path of the reference assembly [./reference.fa]

	Path for the read files
	--samples: path to the directory containing read files [./rawdata]
	--list: the file that lists of the read files [./input.txt]
		less ./input.txt
		test1.1.fq	test1.2.fq
		test2.1.fq	test2.2.fq

	General parameters
	-t: the number of CPU [4]
	-T: the optimum Tm [60]
	-p: the prefix of output file [MultiplexSSR]
	-l: the minimum length of repeat unit [3]
	-n: the minimum number of genotyped individuals [5]
	-a: the minimum number of allele [5]
	-d: the minimum depth of genotype [1]
	-mi: the minimum length of the amplicon [80]
	-ma: the maximum length of the amplicon [470]
	-g: the minimum space between adjacent SSRs [10]
	
	Skip the precedures of genotyping
	-s: skip the procedures of SNP and SSR genotyping, the input file should include SSR.vcf, SNP.vcf and reference.fa [F]

EOQ
}


