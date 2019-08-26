#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: get.all.syn.coord.pl
#
#        USAGE: ./get.all.syn.coord.pl  
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
#      CREATED: 02/14/2019 01:34:57 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;


use Bio::Seq;
use Bio::SeqIO;


my ($allSyn,$wgaDir,$out)=@ARGV;

my @accs=("An-1","C24","Cvi","Eri","Kyo","Ler","Sha");
my %syn;
foreach my $acc (@accs) {
   print "get aln $acc\n";
   my $alnFile = "$wgaDir/$acc/$acc.aligns";
   my $aln = getAln($alnFile);
   my $pairSyn = "$wgaDir/$acc/$acc.wga.syn.block.txt";
   my $syn1 = getPairSYN($allSyn,$pairSyn,$aln);
   $syn{$acc}=$syn1;
   #last;
}

open IN,$allSyn;
open OUT,">$out";
while (<IN>) {
   chomp;
   my @t = (split /\t/);
   print OUT "$t[0]\t$t[1]\t$t[2]";
   #last if ($t[0] eq "Chr2");
   foreach my $acc (@accs) {
   	  print OUT "\t$acc";   	  
   	  if ($syn{$acc}{"$t[0]\t$t[1]\t$t[2]"}{"s"}) {
   	  	my $s = $syn{$acc}{"$t[0]\t$t[1]\t$t[2]"}{"s"};
   	  	print OUT "\t$s";
   	  }else {
   	  	print "##no start position $_\t $acc\n";
   	  	print OUT "\tNA";
   	  }
   	  if ($syn{$acc}{"$t[0]\t$t[1]\t$t[2]"}{"e"}) {
   	  	my $e = $syn{$acc}{"$t[0]\t$t[1]\t$t[2]"}{"e"};
   	  	print OUT "\t$e";
   	  }else {
   	  	print "##no end position $_\t $acc\n";
   	  	print OUT "\tNA";
   	  }
   	  #last;
   }
   print OUT "\n";
}
close OUT;
close IN;



sub getPairSYN {
   my ($allSyn,$pairSyn,$aln)=@_;
   open IN,"intersectBed -a $allSyn -b $pairSyn -wao |";
   print "intersectBed -a $allSyn -b $pairSyn -wao | \n";
   my %syn;
   my %aln=%{$aln};
   while (<IN>) {
   	  chomp;
   	  my @t = (split /\t/);
   	  #last if ($t[0] eq "Chr2");
   	  if ($t[1] >= $t[13] and $t[1] <= $t[14]) {
   	  	my $pairStart = $aln{"$t[12]\t$t[13]\t$t[14]"}{$t[1]};
   	  	if (!$pairStart) {
   	  	   print "#no pair start: $_\n";
   	  	}
   	    if ($syn{"$t[0]\t$t[1]\t$t[2]"}{"s"}  and  $syn{"$t[0]\t$t[1]\t$t[2]"}{"s"} > $pairStart) {
   	    	print "dupStart: $_\tdup\t$pairStart\n";
   	        
   	    }else {
   	       $syn{"$t[0]\t$t[1]\t$t[2]"}{"s"}=$pairStart;	
   	    }   	  	   	  	
   	  }
   	  if ($t[2] >= $t[13] and $t[2] <= $t[14]) {
   	  	my $pairEnd = $aln{"$t[12]\t$t[13]\t$t[14]"}{$t[2]};   	  	
   	  	if (!$pairEnd) {
   	  	  print "no pair end $_\n";
   	  	}
   	  	if ($syn{"$t[0]\t$t[1]\t$t[2]"}{"e"} and $syn{"$t[0]\t$t[1]\t$t[2]"}{"e"} < $pairEnd) {
   	  		print "dupEnd: $_\tdup\t$pairEnd\n";   	  	   
   	  	}else {
   	  		$syn{"$t[0]\t$t[1]\t$t[2]"}{"e"}=$pairEnd;	
   	  	}   	  	  	  	   	  	   	  	   	  	
   	  }
   }     
   close IN;
   foreach my $k (sort keys %syn) {
   	 if ($syn{$k}{"s"} and !$syn{$k}{"e"}) {
   	 	print "#no end: $k\n";
   	 }
   	 elsif (!$syn{$k}{"s"} and $syn{$k}{"e"}) {
   	 	print "#no start: $k\n";
   	 }
   }
   print scalar keys %syn,"\n";
   open IN,$allSyn;
   while (<IN>) {
   	  chomp;
   	  my @t = (split /\t/);
   	  #last if ($t[0] eq "Chr2");
   	  my $k = "$t[0]\t$t[1]\t$t[2]";
   	  print "No start and end \t $k " if (!$syn{$k});
   }
   close IN;
   
   return \%syn; 
}



sub getAln {
   my ($alnFile) = @_;
   open IN,$alnFile;
   my ($chr,$st1,$end1,$st2,$end2,$direc);
   my %aln;
   my @num1; my @num2;
   print "#alnFile\n";
   while (<IN>) {
   	  chomp;
   	  if (/Align/) {
   	  	($chr)=(split /\s+/,$_)[3];
   	  	print "$chr\n";
   	  	#last if ($chr eq "Chr2") ;
   	  }
   	  elsif (/BEGIN/) {
   	  	($st1,$end1,$direc,$st2,$end2) = (split /\s+/ ,$_)[5,7,9,10,12];
   	  	#print  "$st1,$end1,$st2,$end2,$direc \n";   	  	
   	  	@num1=();@num2=();  	  	
   	  }elsif (/^\d/) {
   	  	my ($ast1,$seq1)=(split /\s+/, $_);
   	  	$_=<IN>;
   	  	chomp;
   	  	my ($ast2,$seq2)=(split /\s+/, $_);
   	  	   	  	
   	  	foreach my $i (0..length($seq1)-1) {   	  	   
   	  	   if (substr($seq1,$i,1) ne "\.") {
   	  	   	if (!$num1[-1]) {
   	  	   		push @num1,$st1 
   	  	   	}else {
   	  	   	   push @num1,$num1[-1]+1;	
   	  	   	}   	  	   	     	  	   	  
   	  	   }else {
   	  	   	  if (!$num1[-1]) {
   	  	   	  	push @num1,$st1;
   	  	   	  }else {
   	  	   	    push @num1,$num1[-1];	
   	  	   	  }
   	  	   	  
   	  	   }
   	  	   
   	  	   if (substr($seq2,$i,1) ne "\.") {
   	  	   	 if (!$num2[-1]) {
   	  	   	 	push @num2,$st2;
   	  	   	 }else {
   	  	   	    push @num2,$num2[-1]+1*$direc;	
   	  	   	 }
   	  	   	  
   	  	   }else {
   	  	   	  if (!$num2[-1]) {
   	  	   	  	push @num2,$st2 ;
   	  	   	  }else {
   	  	   	    push @num2,$num2[-1];	
   	  	   	  }   	  	   	  
   	  	   }   	 
   	  	   #if ($st1 == 12091751) {
   	  	   	#  print "num\t$chr\t$num1[-1]\t$num2[-1]\n";
   	  	   #}
   	  	   $aln{"$chr\t$st1\t$end1"}{$num1[-1]}=$num2[-1]; 	
   	  	    #print "> $num1[-1]\t$num2[-1]\n";  
   	  	}
   	  	   	  	    	  	
   	  }else {
   	  	next;
   	  }   	  
   }
   close IN;
   return \%aln;
}
