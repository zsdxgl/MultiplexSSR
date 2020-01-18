#!/usr/bin/perl -w
#use strict;
my $line;my @vol;my $i;my $j;my $k;
my $min_repeats=100;my $max_repeats=0;
my %count_motifs;
my $red_rev;
my ($actual_motif,$actual_motif_a,$reverse_motif,$reverse_motif_a);
open IN, "str_bed_file" or die "$!";
while($line=<IN>){
	chomp $line;
	@vol=split(/\t/,$line);
	$reverse_motif = $actual_motif = $actual_motif_a =$vol[14];
	$reverse_motif =~ tr/ACGT/TGCA/;
	$reverse_motif = reverse $reverse_motif;
	$reverse_motif_a = $reverse_motif;
	for ($j = 0; $j < length ($actual_motif); $j++){
                $actual_motif =~ s/(.)(.*)/$2$1/;
                $reverse_motif =~ s/(.)(.*)/$2$1/;
                $actual_motif_a = $actual_motif if ($actual_motif lt $actual_motif_a);
                $reverse_motif_a = $reverse_motif if ($reverse_motif lt $reverse_motif_a)
        };
        if ($actual_motif_a lt $reverse_motif_a) {$red_rev = "$actual_motif_a/$reverse_motif_a"}
        else {$red_rev = "$reverse_motif_a/$actual_motif_a"}
	$count_motifs{$red_rev}{sprintf("%0.f", $vol[4])}++;
	#print STDERR $count_motifs{$red_rev}{sprintf("%0.f", $vol[4])},"\n";
	if($min_repeats>sprintf("%0.f", $vol[4])){$min_repeats=sprintf("%0.f", $vol[4])};
	if($max_repeats<sprintf("%0.f", $vol[4])){$max_repeats=sprintf("%0.f", $vol[4])};
}
close IN;


my $count=0;
print  "Frequency of classified repeat types (considering sequence complementary)\n\n\nRepeats";
for ($i = $min_repeats; $i <= $max_repeats; $i++) {print  "\t$i"};
print  "\ttotal\n";

foreach $i (sort {length($a)<=>length($b)} keys %count_motifs )
{
	print $i;
	$count=0;
	foreach ($j = $min_repeats; $j <= $max_repeats; $j++)
	{
		if ($count_motifs{$i}{$j}) 
		{
			print  "\t",$count_motifs{$i}{$j};
			$count+=$count_motifs{$i}{$j};
		}
		else {print "\t0"}
	}
	print "\t",$count,"\n";
}
