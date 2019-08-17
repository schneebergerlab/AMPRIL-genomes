#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: calculate.syn.diversity.pl
#
#        USAGE: ./calculate.syn.diversity.pl  
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
#      CREATED: 02/18/2019 12:27:42 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;


my ($allSyn,$alnDir,$lenFileDir,$out)=@ARGV;


#my @accs=("An-1","C24","Col","Cvi","Eri","Kyo","Ler","Sha");
my @accs=("An-1","C24","Col");
## get chrom length
my %len;
foreach my $acc (@accs) {
  my $in = "$lenFileDir/$acc.leng.txt";
  $len{$acc} = getLeng($in);
}

## get coordinates of syteny-all regions
my $ref = "Col";
open IN,$allSyn; ## 1 based coord,
my %coord;
my $ch;
while (<IN>) {
  chomp;
  my @t = (split /\t/);  
  if (!$ch) {
  	my @ts=($t[0],-1,0);
  	foreach my $i (0..$#accs) {
  	  push @ts, $accs[$i];
  	  push @ts , -1;
  	  push @ts, 0;  	  
  	}
  	push @{$coord{$t[0]}},[@ts]; #incorporate a start "SYN"
  	push @{$coord{$t[0]}},[@t];
  	$ch = $t[0];
  }elsif ($ch ne $t[0]) {
  	my @te= ($ch, $len{$ref}{$ch}+1, $len{$ref}{$ch}+1);  	
  	my @ts=($t[0], -1, 0);
  	foreach my $i (0..$#accs) {
  	  @ts = (@ts, $accs[$i], -1 ,0) ;
  	  @te = (@te, $accs[$i], $len{$accs[$i]}{$ch}+1, $len{$accs[$i]}{$ch}+1); #incorporate a end "SYN"	    	  
  	}
  	push @{$coord{$ch}},[@te];
  	$ch = $t[0];
  	push @{$coord{$ch}},[@ts]; 
  	push @{$coord{$ch}},[@t];;
  }else {
    push @{$coord{$t[0]}},[@t];	
  }  
}
my @te= ($ch,$len{$ref}{$ch}+1, $len{$ref}{$ch}+1);  
foreach my $i (0..$#accs) {  	 
  @te = (@te, $accs[$i], $len{$accs[$i]}{$ch}+1, $len{$accs[$i]}{$ch}+1);	    	  
}
push @{$coord{$ch}},[@te];  	
close IN;

##### calculate "SYN diversity score for each site"
my %aln; my $rIx=2;
my %alnRef; ## hash the alignment positions of ref and acc in syn aln blocks 
my %trans; ## used to transfer the coordinate
my %div; ##synteny diversity
my $comNum = 0;
print "##start to hash the  SV-type for each postion of each  pairwise alignment and the Col-ref based pairwise alignment \n";
foreach my $i (0..$#accs-1) {
  foreach my $j ($i+1..$#accs) {  	
  	 my $acc1 = $accs[$i];
  	 my $acc2 = $accs[$j];
  	 print "# $acc1\t$acc2\n";
  	 $comNum+=1;  	 
  	 
  	 if ($acc1 eq $ref) {  	 	  	
  	 	my $alnFile = "$alnDir/$acc1/$acc2/$acc2.aligns";
  	    my $blkFile = "$alnDir/$acc1/$acc2/$acc2.wga.block.txt";  	  	  	 	
  	 	my $blk = getBlk($blkFile);
  	 	my ($alnRef,$trans) = getAlnRef($alnFile,$blk);
  	 	$alnRef{$acc1}{$acc2} = $alnRef;
  	 	$trans{$acc1}{$acc2} = $trans;
  	 }elsif ($acc2 eq $ref) {
  	 	my $alnFile = "$alnDir/$acc2/$acc1/$acc1.aligns";
  	    my $blkFile = "$alnDir/$acc2/$acc1/$acc1.wga.block.txt";  	  	  	 	
  	 	my $blk = getBlk($blkFile);
  	 	my ($alnRef,$trans) = getAlnRef($alnFile,$blk);
  	 	$alnRef{$acc2}{$acc1} = $alnRef;
  	 	$trans{$acc2}{$acc1} = $trans;
  	 }else {
  	 	my $alnFile = "$alnDir/$acc1/$acc2/$acc2.aligns";
  	    my $blkFile = "$alnDir/$acc1/$acc2/$acc2.wga.block.txt";  	  	  	 	
  	 	my $blk = getBlk($blkFile);
  	 	my ($aln1,$aln2) = getAln($alnFile,$blk);
  	 	$aln{$acc1}{$acc2} = $aln1;
  	 	$aln{$acc2}{$acc1} = $aln2;  	 	
  	 }
  }  
}
print "#combinations $comNum\n";

open OUT,">$out";

foreach my $chr (sort keys %coord) {  
  my @syns = @{$coord{$chr}};       
  print "##SYN $chr\t$#syns\n";
  foreach my $i (1..$#syns) {
   my ($refSt,$refEnd)=($syns[$i-1][3+$rIx*3+2] + 1, $syns[$i][3 + $rIx*3 + 1] - 1); # 1 based coord,
   foreach my $pos ($refSt..$refEnd) {   	       	  
   	 $div{$chr}{$pos} = [(0,0)]  ## syn vs non-syn
   }
   foreach my $j (0..$#accs-1) {   	   	
   	 foreach my $k ($j+1..$#accs) {
   	    my $acc1 = $accs[$j];
   	    my $acc2 = $accs[$k];
   	    my ($start1,$end1) = ($syns[$i-1][3+$j*3+2] + 1, $syns[$i][3+$j*3+1]-1);
   	    my ($start2,$end2) = ($syns[$i-1][3+$k*3+2] + 1, $syns[$i][3+$k*3+1]-1);
   	    if ($acc1 eq $ref) {
   	      foreach my $pos ($start1..$end1) {   	      	
   	      	if ($alnRef{$acc1}{$acc2}{$chr}{$pos}) {
   	      		$div{$chr}{$pos}[0]+= 1;
   	      	}else {
   	      		$div{$chr}{$pos}[1]+= 1;
   	      	}
   	      }
   	    }elsif ($acc2 eq $ref) {
   	      foreach my $pos ($start2..$end2) {   	      	
   	      	if ($alnRef{$acc2}{$acc1}{$chr}{$pos}) {
   	      		$div{$chr}{$pos}[0]+= 1;
   	      	}else {
   	      		$div{$chr}{$pos}[1]+= 1;
   	      	}
   	      }   	    	
   	    }else {
   	      
   	       foreach my $pos ($refSt..$refEnd) {
   	       	  if (! $trans{$ref}{$acc1}{$chr}{$pos} and ! $trans{$ref}{$acc2}{$chr}{$pos}) { ## both del
   	       	  	 $div{$chr}{$pos}[0] += 1
   	       	  }elsif (!$trans{$ref}{$acc1}{$chr}{$pos} ) {
   	       	  	 $div{$chr}{$pos}[1] += 1;
   	       	  }elsif (! $trans{$ref}{$acc2}{$chr}{$pos}) {
   	       	  	 $div{$chr}{$pos}[1] += 1;
   	       	  }else {
   	       	  	my @pos1 = sort {$a<=>$b} keys %{$trans{$ref}{$acc1}{$chr}{$pos}};
   	       	  	my @pos2 = sort {$a<=>$b} keys %{$trans{$ref}{$acc2}{$chr}{$pos}};
   	       	  	my $syn = 0;
   	       	  	#print "aa: $#pos1\t$#pos2\t$pos1[0]\t$pos2[0]\n";
   	       	  	foreach my $pos1 (@pos1) {
   	       	  		#next if ($pos1 < $start1 or $pos1 > $end1);
   	       	  		foreach my $pos2 (@pos2) {
   	       	  		   #next if ($pos2 < $start2 or $pos2 > $end2);
   	       	  		   if ($aln{$acc1}{$acc2}{$chr}{$pos1} and $aln{$acc2}{$acc1}{$chr}{$pos2}  ) {
   	       	  		   	  $syn = 1;
   	       	  		   	  #print "#yy $acc1\t$acc2\t$chr\t$pos\t$pos1\t$pos2\n";
   	       	  		   }	
   	       	  		}
   	       	  	}
   	       	  	if ($syn == 1) {
   	       	  	   $div{$chr}{$pos}[0] += 1;
   	       	  	}else {
   	       	  	   $div{$chr}{$pos}[1] += 1;
   	       	  	}
   	       	  }   	       	  
   	       }
   	    }   	    	  	 
     } 	   	  	 
   }
   foreach my $pos ($refSt..$refEnd) {   	       	     	 
   	  my $t = $div{$chr}{$pos}[0]+$div{$chr}{$pos}[1];
   	  print "$chr\t$pos\t$t\t$div{$chr}{$pos}[0]\t$div{$chr}{$pos}[1]\n" if ( $t != $comNum );
   	  my $p = sprintf("%.3f", $div{$chr}{$pos}[1] / $t);
   	  print OUT "$chr\t$pos\t$t\t$p\n";
   }  
   ($refSt,$refEnd)=($syns[$i][3+$rIx*3+1], $syns[$i][3 + $rIx*3 + 2]);
   foreach my $pos ($refSt..$refEnd) {   	       	     	
   	  print OUT "$chr\t$pos\t$comNum\t1.000\n";
   }       
 }     
 #last; ##TEST
}
close OUT;




sub getAln {
   my ($alnFile,$blk) = @_;
   open IN,$alnFile;
   my ($chr,$st1,$end1,$st2,$end2,$direc,$sv);
   my %aln1; my %aln2;
   my %blk = %{$blk};
   my @num1; my @num2;
   
   print "#alnFile\n";
   while (<IN>) {
   	  chomp;
   	  if (/Align/) {
   	  	($chr)=(split /\s+/,$_)[3];
   	  	$chr=~s/chr/Chr/g;
   	  	print "$chr\n";
   	  	#last if ($chr eq "Chr2") ; # TEST
   	  }
   	  elsif (/BEGIN/) {
   	  	($st1,$end1,$direc,$st2,$end2) = (split /\s+/ ,$_)[5,7,9,10,12];
   	  	#print  "$st1,$end1,$st2,$end2,$direc \n";   	  	
   	  	if ($blk{$chr}{"$st1\t$end1\t$st2\t$end2"}) {
   	  	  $sv = $blk{$chr}{"$st1\t$end1\t$st2\t$end2"};
   	  	}else {
   	  	   $sv = "NA";
   	  	}
   	  	@num1=();@num2=();  	  	
   	  }elsif (/^\d/) {
   	  	next if ($sv ne "SYN");
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
   	  	   #$aln{"$chr\t$st1\t$end1"}{$num1[-1]}=$num2[-1]; 	
   	  	    #print "> $num1[-1]\t$num2[-1]\n";
   	  	    
   	  	   $aln1{$chr}{$num1[-1]} = $sv;
   	  	   $aln2{$chr}{$num2[-1]} = $sv;   
   	  	}
   	  	   	  	    	  	
   	  }else {
   	  	next;
   	  }   	  
   }
   close IN;
   return (\%aln1,\%aln2);
}

sub getBlk {
  my ($in) = @_;
  open IN,$in;
  my %blk;
  while (<IN>) {
  	chomp;
  	my @t = (split /\t/);
  	next if ($t[0] ne $t[3]);
  	$blk{$t[0]}{"$t[1]\t$t[2]\t$t[4]\t$t[5]"}=$t[-1];  	
  }
  return \%blk;
}

sub getAlnRef {
   my ($alnFile,$blk) = @_;
   open IN,$alnFile;
   my ($chr,$st1,$end1,$st2,$end2,$direc,$sv);
   my %aln1; my %aln2;
   my %blk = %{$blk};
   my @num1; my @num2;
   
   print "#alnFile\n";
   while (<IN>) {
   	  chomp;
   	  if (/Align/) {
   	  	($chr)=(split /\s+/,$_)[3];
   	  	print "$chr\n";
   	  	#last if ($chr eq "Chr2") ; #TEST
   	  }
   	  elsif (/BEGIN/) {
   	  	($st1,$end1,$direc,$st2,$end2) = (split /\s+/ ,$_)[5,7,9,10,12];
   	  	#print  "$st1,$end1,$st2,$end2,$direc \n";   	  	
   	  	if ($blk{$chr}{"$st1\t$end1\t$st2\t$end2"}) {
   	  	  $sv = $blk{$chr}{"$st1\t$end1\t$st2\t$end2"};
   	  	}else {
   	  	   $sv = "NA";
   	  	}
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
   	  	   #$aln{"$chr\t$st1\t$end1"}{$num1[-1]}=$num2[-1]; 	
   	  	    #print "> $num1[-1]\t$num2[-1]\n";
   	  	    
   	  	   #$aln{$acc1}{$acc2}{$chr}{$num1[-1]} = $sv;
   	  	   #$aln{$acc2}{$acc1}{$chr}{$num2[-1]} = $sv;   
   	  	   #$aln{$acc1}{$acc2}{$chr}{$num1[-1]}{$num2[-1]} = $sv
   	  	   $aln2{$chr}{$num1[-1]}{$num2[-1]} = $sv; ## for coordinate translation
   	  	   $aln1{$chr}{$num1[-1]}= $sv if ($sv eq "SYN");
   	  	}
   	  	   	  	    	  	
   	  }else {
   	  	next;
   	  }   	  
   }
   close IN;
   return (\%aln1,\%aln2);
}


sub getLeng {
  my ($in) = @_;
  open IN,$in;
  my %len;
  while (<IN>) {
  	 chomp;
  	 my @t = (split /\t/);
  	 $len{$t[0]}=$t[1];
  }
  close IN;
  return \%len;
}
