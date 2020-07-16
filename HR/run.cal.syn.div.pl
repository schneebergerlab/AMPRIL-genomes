#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: run.cal.syn.div.pl
#
#        USAGE: ./run.cal.syn.div.pl  
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
#      CREATED: 02/19/2019 11:25:34 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my ($allSyn,$alnDir,$lenFileDir,$script,$outdir) = @ARGV;

my @accs=("An-1","C24","Col","Cvi","Eri","Kyo","Ler","Sha");
my $n=5;
my $chr="Chr";
my $num=0;
for my $i (0..$#accs-1) {
  for my $j ($i+1..$#accs) {
  	my $acc1 = $accs[$i];
  	my $acc2 = $accs[$j];
  	for my $k (1..$n) {
  	  my $out = "$outdir/$acc1.$acc2.Chr$k.pairwise.syn.txt";
  	  my $log = "$outdir/np.$acc1.$acc2.Chr$k.log";
  	  my $chrom = $chr.$k;
  	  print "perl $script $allSyn $alnDir $lenFileDir $acc1 $acc2 $chrom $out > $log \n" ;
      system("perl $script $allSyn $alnDir $lenFileDir $acc1 $acc2 $chrom $out > $log &" );
      $num+=1;
      if ($num==10) {
      	sleep(800);
      	$num=0;
      }  		
  	}
  }
} 

