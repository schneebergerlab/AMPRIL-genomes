#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: syn.div.merge.pl
#
#        USAGE: ./syn.div.merge.pl  
#
#  DESCRIPTION: a
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 02/19/2019 02:38:46 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my ($indir,$chr,$lengFile,$out)=@ARGV;

my @accs=("An-1","C24","Col","Cvi","Eri","Kyo","Ler","Sha");

open IN,$lengFile;
my @div;
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  next if ($t[0] ne $chr);
  @div = (0) x $t[1];
}
close IN;

my $comb=0;
for my $i (0..$#accs-1) {
  for my $j ($i+1..$#accs) {
  	my $acc1 = $accs[$i];
  	my $acc2 = $accs[$j];
  	$comb+=1;
  	my $file = "$indir/$acc1.$acc2.$chr.pairwise.syn.txt";
  	open IN,$file;
  	while (<IN>) {
  	  chomp;
  	  my @t = (split /\t/);
  	  $div[$t[1]-1]+=1 if $t[3]==1;
  	}
  	close IN;
  }
} 
open OUT,">$out";
foreach my $i (0..$#div) {
  my $pos = $i+1;
  my $p = sprintf("%.3f", $div[$i]/$comb);
  print OUT "$chr\t$pos\t$comb\t$div[$i]\t$p\n";
}
close OUT;

