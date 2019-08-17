#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: calculate.syn.diversity.window.pl
#
#        USAGE: ./calculate.syn.diversity.window.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 02/18/2019 06:29:48 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;


my ($inFile,$win,$step,$out) = @ARGV;


open IN,$inFile;
my %div; my $leng;
my $chr; my $ncom;
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  my $w = int(($t[1]-1)/$step);
  $div{$w}+=$t[3];
  $leng = $t[1];
  $chr = $t[0];  
  $ncom = $t[2];
}
close IN;

my $start=1;
my $end=$start+$win-1;
my @wdiv;
my $ns = $win/$step;

open OUT,">$out";
my $div = 0;
my $n=0;
while ($end <= $leng) {  
  my $ss = int($start/$step);
  $div = 0;
  
  foreach my $i (0..$ns-1) {
  	 $div += $div{$i+$ss};
  } 
  #push @wdiv,$div;
  my $divPer = sprintf( "%.3f",$div/($win*$ncom) );
  print OUT "$chr\t$start\t$end\t$n\t$divPer\n";
  $start += $step;
  $end = $start + $win - 1;
  $n++;
}

my $divPer = sprintf("%.3f",$div/( ( $leng-$start+1 ) * $ncom ) );
print OUT "$chr\t$start\t$leng\t$n\t$divPer\n";


