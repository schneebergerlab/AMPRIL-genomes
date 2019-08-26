#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: R.gene.cluster.wga.ortho.scheme.pl
#
#        USAGE: ./R.gene.cluster.wga.ortho.scheme.pl  
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
#      CREATED: 03/15/2019 04:18:35 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use SVG;
use FileHandle;
use Bio::Seq;
use Bio::SeqIO;
#use Pod::Usage;
#use Getopt::Long;

my ($allSynFile, $rCluBed, $alnDir, $ortFile, $geneBedDir, $annDir, $win , $outDir)=@ARGV;

system("mkdir -p $outDir");
#my $win = 10000;
### step 1: get the alignments in-between syn-all regions
print "windowBed -w $win -b $allSynFile -a $rCluBed |sort -k6,6 -k7,7 -k8,8n -k9,9n | \n";
open IN, "windowBed -w $win -b $allSynFile -a $rCluBed |sort -k6,6 -k7,7 -k8,8n -k9,9n |";
my %clu; 
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  push @{$clu{$t[5]}},[@t];
}
close IN;
my @accs2=("An-1","C24","Cvi","Eri","Kyo","Ler","Sha");
my @accs1 =("An-1","C24","Col","Cvi","Eri","Kyo","Ler","Sha");
my @accs3=("Col","An-1","C24","Cvi","Eri","Kyo","Ler","Sha");
my %reg;

#define the boundary of all-syn regions
my %beds;
foreach my $acc (@accs1) {
  #my $fh = FileHandle->new("> $outDir/$acc.bed");
  my $fh = FileHandle->new();
  open ($fh, "|sort -k1,1 -k2,2n -k3,3n > $outDir/$acc.bed");
  #$FH{"filtered"} = $fh;
  $beds{$acc} = $fh;
  #print $fh,"\n";
}

my %max;
my $s=9;
foreach my $clu (sort keys %clu) {
  my @tmp= @{$clu{$clu}};
  my @ts = @{$tmp[0]};
  my @te = @{$tmp[-1]};
  my $chr = $ts[0];
  my $maxA=0;
  foreach my $i (0..$#accs1) {
  	 my $acc = $accs1[$i];
  	 $reg{$clu}{$acc} = [($chr, $ts[$i*3+$s+1], $te[$i*3+$s+2])];  	 
  	 if ( $te[$i*3+$s+2] - $ts[$i*3+$s+1] > $maxA ) {
  	 	$maxA = $te[$i*3+$s+2] - $ts[$i*3+$s+1];
  	 }
  	 my $fhp = FileHandle->new();
  	 $fhp = $beds{$acc};
  	 print $fhp "$chr\t$ts[$i*3+$s+1]\t$te[$i*3+$s+2]\t$clu\n";
  }
  $max{$clu} = $maxA;
}

foreach my $acc (@accs1) {
  my $fhp = FileHandle->new();
  $fhp = $beds{$acc};
  close $fhp;
}
 	
#get the alignments
my %aln;  	
foreach my $clu (sort keys %reg) {
  system("mkdir -p $outDir/$clu");
  my $aln = getAln($reg{$clu}, $clu, \@accs3, "$outDir/$clu");
  $aln{$clu}=$aln;
}

#exit;

## step 2: get orthologs and R genes in the regions
my $ort = getOrtho($ortFile,\@accs2);
#my %ort = %{$ort};

my %genesTmp; my %rgenes; my %allrgenes;
foreach my $i (0..$#accs1) {
  my $acc = $accs1[$i];    
  my $seq = getProt("$geneBedDir/$acc.prot.v2.5.fasta");
  my $genes = getGenes("$geneBedDir/$acc.protein-coding.genes.v2.5.bed","$outDir/$acc.bed", $seq,$outDir,$acc);
  $genesTmp{$acc}=$genes;
  
  #
  my $geneAnn = "$annDir/$acc.version2.5.RGA.bed";
  my $ann = getRgenes($geneAnn);  
  $rgenes{$acc} = $ann;
  %allrgenes = (%allrgenes,%{$ann});
}
my %genes;
foreach my $acc (sort keys %genesTmp) {
   foreach my $clu (sort keys %{$genesTmp{$acc}}) {
   	  $genes{$clu}{$acc} = $genesTmp{$acc}{$clu};
   }
}


foreach my $clu (sort keys %clu) {
   
   print $clu,"\n";
   if (!$genes{$clu}) {
   	  print "#no genes in $clu \n";
   	  next;
   }
   
   my $ortRegFile = "$outDir/$clu/prot/Results_Jul30/Orthogroups.txt";
   if (! -e $ortRegFile) {
   	  print "$ortRegFile does not exist\n";
   	  system("orthofinder -t 10 -a 2 -f $outDir/$clu/prot");   	   
   }
   my $outFile = "$outDir/$clu/local.orthogroups.txt";
   grToCSV($ortRegFile,\@accs1,$outFile);

   my $ort = getOrtho($outFile,\@accs2);
   my ($gr,$rgr) = getOrtho2($outFile,\%allrgenes);

   my $outSVG = "$outDir/$clu/$clu.AMPRIL.svg";
   my $outPDF = "$outDir/$clu/$clu.AMPRIL.pdf";
      
   drawFig2($reg{$clu}, $aln{$clu}, $max{$clu}, $genes{$clu}, $ort, \%rgenes, \@accs3, $outSVG,0, $gr,$rgr);   
   print "cairosvg -f pdf -o $outPDF $outSVG\n";
   system("cairosvg -f pdf -o $outPDF $outSVG");
   
   $outSVG = "$outDir/$clu/$clu.AMPRIL.ortR.svg";
   $outPDF = "$outDir/$clu/$clu.AMPRIL.ortR.pdf";
   drawFig2($reg{$clu}, $aln{$clu}, $max{$clu}, $genes{$clu}, $ort, \%rgenes, \@accs3, $outSVG,1, $gr,$rgr);   
   print "cairosvg -f pdf -o $outPDF $outSVG\n";
   system("cairosvg -f pdf -o $outPDF $outSVG");
   
   $outSVG = "$outDir/$clu/$clu.AMPRIL.aln.svg";
   $outPDF = "$outDir/$clu/$clu.AMPRIL.aln.pdf";
   drawFig2($reg{$clu}, $aln{$clu}, $max{$clu}, $genes{$clu}, $ort, \%rgenes, \@accs3, $outSVG,2, $gr,$rgr);   
   print "cairosvg -f pdf -o $outPDF $outSVG\n";
   system("cairosvg -f pdf -o $outPDF $outSVG");
   #exit;	
}

exit;

foreach my $clu (sort keys %clu) {
   my $outSVG = "$outDir/$clu/$clu.AMPRIL.svg";
   my $outPDF = "$outDir/$clu/$clu.AMPRIL.pdf";
      
   drawFig($reg{$clu}, $aln{$clu}, $max{$clu}, $genes{$clu}, $ort, \%rgenes, \@accs3, $outSVG,0);   
   print "cairosvg -f pdf -o $outPDF $outSVG\n";
   system("cairosvg -f pdf -o $outPDF $outSVG");
   
   $outSVG = "$outDir/$clu/$clu.AMPRIL.ortR.svg";
   $outPDF = "$outDir/$clu/$clu.AMPRIL.ortR.pdf";
   drawFig($reg{$clu}, $aln{$clu}, $max{$clu}, $genes{$clu}, $ort, \%rgenes, \@accs3, $outSVG,1);   
   print "cairosvg -f pdf -o $outPDF $outSVG\n";
   system("cairosvg -f pdf -o $outPDF $outSVG");
   
   $outSVG = "$outDir/$clu/$clu.AMPRIL.aln.svg";
   $outPDF = "$outDir/$clu/$clu.AMPRIL.aln.pdf";
   drawFig($reg{$clu}, $aln{$clu}, $max{$clu}, $genes{$clu}, $ort, \%rgenes, \@accs3, $outSVG,2);   
   print "cairosvg -f pdf -o $outPDF $outSVG\n";
   system("cairosvg -f pdf -o $outPDF $outSVG");
   #exit;	
}


sub grToCSV { ## change the output of orthoFinder into tab file including each accession infor
  my ($in,$accs,$out)=@_;
  open IN,$in;
  my @accs1 = @{$accs};
 
  open OUT,">$out";
  while (<IN>) {
  	 chomp;
  	 my @t = (split / /,$_);
  	 my $gr = $t[0];
  	 shift @t;
  	 my %gene;
  	 #print "$in\t$gr\t$t[0]\n";
  	 foreach my $gene (@t) {
  	 	if ($gene=~/AT\d/) {
  	 	  push @{$gene{"Col"}},$gene;	
  	 	}elsif ($gene=~/ATAN/) {
  	 	  push @{$gene{"An-1"}},$gene;	
  	 	}
  	 	elsif ($gene=~/ATC24/) {
  	 	  push @{$gene{"C24"}},$gene;	
  	 	}
  	 	elsif ($gene=~/ATCVI/) {
  	 	  push @{$gene{"Cvi"}},$gene;	
  	 	}
  	 	elsif ($gene=~/ATERI/) {
  	 	  push @{$gene{"Eri"}},$gene;	
  	 	}
  	 	elsif ($gene=~/ATKYO/) {
  	 	  push @{$gene{"Kyo"}},$gene;	
  	 	}
  	 	elsif ($gene=~/ATLER/) {
  	 	  push @{$gene{"Ler"}},$gene;	
  	 	}
  	 	elsif ($gene=~/ATSHA/) {
  	 	  push @{$gene{"Sha"}},$gene;	
  	 	}
  	 }
  	 my @line; push @line,$gr;
  	 foreach my $acc ( @accs1) {
  	 	if ($gene{$acc}) {
  	 	  push @line, join(",",@{$gene{$acc}});  	 	  
  	 	}else {
  	 		push @line,"-----";
  	 	}  	 	
  	 }
  	 my $line = join("\t",@line);
  	 print OUT "$line\n";
  }
  close IN;close OUT;
	
}


sub getGenes {
	my ($geneBed, $regBed,$seq,$outdir,$acc)=@_;
	my %genes; 
	open IN, "intersectBed -a $regBed -b $geneBed -wao |";
	print "intersectBed -a $regBed -b $geneBed -wao | \n";
	my %regProt; my %seq=%{$seq};
	while (<IN>) {
		chomp;
		my @t = (split /\t/);
		next if ($t[-1]==0);
		next if ($t[-1] < $t[-4]-$t[5]-100);
		#print "$t[2]\t$t[-2]\n";exit;
		my $id= (split /\;/,$t[-2])[0];
		$id = (split /\=/,$id)[1];
		#push @genes, [($id,$t[-5],$t[-4],$t[-3])];
		$genes{$t[3]}{$id}=[($id,$t[-5],$t[-4],$t[-3])];
		$regProt{$t[3]}{$id}=$seq{$id};
	}
	close IN;
	#my $n = scalar keys %genes;
	#print "total genes $geneBed in the region: $n \n";
	foreach my $reg (sort keys %regProt) {
	  system("mkdir -p $outdir/$reg/prot");
	  open OUT,">$outdir/$reg/prot/$acc.fasta";
	  foreach my $id (sort keys %{$regProt{$reg}}) {
	  	print OUT ">$id\n$regProt{$reg}{$id}\n";
	  }
	  close OUT;
	}
	
	return \%genes;
}

sub getProt {
	my ($in)=@_;
	my %seq;
	my $seqin=new Bio::SeqIO(-file=>$in,-format=>"fasta");
	while (my $seqobj=$seqin->next_seq()) {
		my $id=$seqobj->id();
		($id)=(split /\./,$id)[0];
		my $seq=$seqobj->seq();
		$seq{$id}=$seq;
	}
	return \%seq;
}


sub getRgenes {
	my ($geneAnn)=@_;
	open IN, "$geneAnn";
	my %genes;
	while (<IN>) {
		chomp;
		my @t = (split /\t/);
		#$t[0] = (split /\./,$t[0])[0];
		#next if ($t[5] ne "NBS-LRR");
		if ($t[5] eq "NBS-LRR") {
		  if ($t[6] eq "TNL" or $t[6] eq "TN") {
			$genes{$t[4]} = 1;  
		  }elsif ($t[6] eq "CNL" or $t[6] eq "CN") {
			$genes{$t[4]} = 2;			  
		  }else {
			$genes{$t[4]} = 3;
		  }
		}elsif ($t[5] eq "RLK") {
			$genes{$t[4]} = 4
		}elsif ($t[5] eq "RLP") {
			$genes{$t[4]} = 5;
		}else {
			$genes{$t[4]} = 6;
		}
		
	}
	close IN;
	return \%genes;
}


sub getAln {
	my ($reg,$clu, $accs,$outDir) = @_;
	my %reg = %{$reg};
	my %aln;
	my @accs = @{$accs};
	foreach my $i (0..$#accs-1) {
   	  my $acc1 = $accs[$i];
   	  my $acc2 = $accs[$i+1];
   	  my $alnFile = "$alnDir/$acc1/$acc2/$acc2.wga.block.txt";
   	  open IN,$alnFile;
   	  open OUT,">$outDir/$clu.$acc1.$acc2.wga.aln.txt";
      my $n = 1;
      my $chr = $reg{$acc1}[0];      
      while (<IN>) {
   	    chomp;
   	    my @t = (split /\t/);
   	    next if ($t[0] ne $chr);
   	    next if ($t[-1]=~/DUP/);
   	    next if ($t[-1] eq "CTX");
   	    next if ($t[-1] eq "ITX" and $t[5] < $reg{$acc2}[1]);
   	    next if ($t[-1] eq "ITX" and $t[4] > $reg{$acc2}[2]);
   	   
   	    if ($reg{$acc1}[1] >= $t[1] and $reg{$acc1}[1] <= $t[2]) {
   	  	  push @{$aln{$acc1}{$acc2}},[@t];
   	  	  print OUT "$_\n";
   	    }elsif ($t[1]>=$reg{$acc1}[1] and $t[2] <= $reg{$acc1}[2]) {
   	  	  push @{$aln{$acc1}{$acc2}},[@t];
   	  	  print OUT "$_\n";
   	    }elsif ($reg{$acc1}[2] >= $t[1] and $reg{$acc1}[2] <= $t[2]) {
   	  	  push @{$aln{$acc1}{$acc2}},[@t];
   	  	  print OUT "$_\n";
   	    }
      }
      close IN;      
      close OUT;
    }
    return \%aln;
}


sub getOrtho {
   my ($infile,$accs) = @_;
   open IN,$infile;
   my @accs = @{$accs};
   my %ort;
   while (<IN>) {
   	  chomp;
   	  my @t = (split /\t/);
   	  shift @t;
   	  foreach my $i (0..$#accs) {
   	  	#next if ($i==2);
   	  	if ($i==0) {
   	  	   next if ($t[2] eq "-----");
   	  	   next if ($t[$i] eq "-----");
   	  	   my @genes = (split /\,/,$t[$i]);   	  	   
   	  	   my @cols = (split /\,/,$t[2]);
   	  	   foreach my $g (@genes) {
   	  	   	  $ort{$g}= [@cols];
   	  	   }
   	  	}else {   	  	   
   	  	  my @genes ; my @prev;   	  	   
   	  	   if ($i == 2) {
   	  	   	  next if ($t[$i+1] eq "-----");
   	  	   	  next if ($t[$i-1] eq "-----");
   	  	   	  @genes = (split /\,/,$t[$i+1]);
   	  	      @prev = (split /\,/,$t[$i-1]);
   	  	   }elsif ($i>2) {
   	  	   	  #print "xxx:$accs[$i]\t$t[$i+1]\t$t[$i]\n";
   	  	   	  next if ($t[$i+1] eq "-----");
   	  	   	  next if ($t[$i] eq "-----");
   	  	   	  @genes = (split /\,/,$t[$i+1]);
   	  	      @prev = (split /\,/,$t[$i]);
   	  	   }else {
   	  	   	  next if ($t[$i] eq "-----");
   	  	   	  next if ($t[$i-1] eq "-----");
   	  	   	  @genes = (split /\,/,$t[$i]);
   	  	      @prev = (split /\,/,$t[$i-1]);
   	  	   }
   	  	   
   	  	   foreach my $g (@genes) {
   	  	   	  $ort{$g}= [@prev];
   	  	   	  
   	  	   }	
   	  	}
   	  }
   }
   #print $ort{"ATAN-4G46470"}, $ort{"ATCVI-4G43710"},"\n";
   #print "test ort:", $ort{"ATAN-5G69390"}[0],"\t", $ort{"ATCVI-5G67340"}[0],"\n";
   return \%ort;   
}

sub getOrtho2 {
	my ($infile,$rgene) = @_;
	my %rgene=%{$rgene};
	open IN,$infile;
	my %rgr;my %gr; my %genes;
	my $k = 0;
	while (<IN>) {
		chomp;
		my @t = (split /\t/);
		my $gr=$t[0];
		shift @t;
		my $flag = 0;
		foreach my $genes (@t) {
			my @genes = (split /\,/,$genes);
			foreach my $gene (@genes) {
			  $gr{$gene}=$gr;			
			  $flag = 1 if $rgene{$gene};	
			}						 
		}
		if ($flag==1) {		  
		  $rgr{$gr}=$k+1;
		  $k+=1;
		}else {
		  $rgr{$gr}=0;	
		}
	}
	close IN;
	return (\%gr, \%rgr);
}


sub drawFig {
	# https://www.w3.org/TR/SVG/paths.html#PathDataLinetoCommands
	# https://metacpan.org/pod/Math::Bezier
	# https://metacpan.org/pod/SVG#path	
   my ($reg, $aln, $maxSize, $genes, $ort, $rgenes , $accs, $outSVG,$flag) = @_;
   
   my @accs = @{$accs};
   my %aln = %{$aln};
   my %reg = %{$reg};
   	
   my %genes = %{$genes};
   my %ort = %{$ort};
   my %rgenes = %{$rgenes};
   
   my $unit=500; #scale	
   my $width = $maxSize/$unit + 120 ;
   my $height = ($#accs+1)*40 + 120;		
   my $svg= SVG->new(width=>$width,height=>$height);
   print "#figure size $outSVG $width  x  $height\n";     
	
	my $gold="rgb(255,215,0)";
	my $red="rgb(255,0,0)";
	my $lred="rgb(128,0,0)";
	my $blue="rgb(0,0,255)";
	my $lblue="rgb(0,0,128)";
	my $green="rgb(0,255,0)";
	my $lgreen="rgb(0,128,0)";
	my $gray="rgb(24,24,24)";
	my $llgray="rgb(196,196,196)";
	my $lgray="rgb(97,96,96)";
	my $black="rgb(0,0,0)";
	my $dgreen="rgb(53,128,0)";
    my $vvlblue = "rgb(239,243,255)";
    my $orange = "rgb(255, 127, 0 )";
    my $puple = "rgb(128,0,128)";
	##top fig for Col
	
	#draw the Col reg 
	my $x0 = 40;
	my $y0 = 40; 	
	my $top = 25; # distance between up and down regions
	my $hei = 10; #rectangle height for gene
	my $h = 2 ; ##  rectangle height for the region
	my $up = int(($hei-$h)/2); ## the height of up region of a gene
	my $fsize = 3 ; # font size  
	
    my $endRef = $reg{"Col"}[2];
	my $startRef = $reg{"Col"}[1];
	$svg->text(id=>"Col-0", x=>$x0-30, y=>$y0, "font-family"=>"Arial","font-size"=>$fsize*2, -cdata => "Col", );
	$svg->rectangle(x=>$x0,y=>$y0,width=>($endRef-$startRef)/$unit, height=>$h,style=>{'fill'=>$lgray,'stroke'=>$lgray, 'stroke-width'=>0,});
		
	#start position
	$svg->text(x=>$x0, y=>$y0-20, "font-family"=>"Arial","font-size"=>$fsize, -cdata => $startRef, style=>{'writing-mode'=>'tb','text-anchor'=>"middle",});
	#end position
	$svg->text(x=>$x0+($endRef-$startRef)/$unit, y=>$y0-20, "font-family"=>"Arial","font-size"=>$fsize, -cdata => $endRef, style=>{'writing-mode'=>'tb','text-anchor'=>"middle",});
	
	
    my %pos;
    
	foreach my $gene (sort keys %{$genes{"Col"}}) {
		my ($id,$start,$end,$direc) = @{$genes{"Col"}{$gene}};
		my $color = $lblue;      
		print "gene $gene\n";  		
		if ($rgenes{"Col"}{$gene}) {					
			if ($rgenes{"Col"}{$gene}==1) {
			  	$color = $red; ##TNL,TN
			}elsif ($rgenes{"Col"}{$gene}==2 ) {
			    $color = $orange; ##CNL,CN
			}elsif ($rgenes{"Col"}{$gene}==3){
			    $color = $green; #other NBS-LRR
			}elsif ($rgenes{"Col"}{$gene}==4){
			    $color = $gold; ##RLK
			}	elsif ($rgenes{"Col"}{$gene}==5){
			    $color = $gray; ##RLP
			}else {
				$color = $puple; ##TM-CC
			}						
			$svg->text(x=>$x0+($start-$startRef)/$unit, y=>$y0-20, "font-family"=>"Arial","font-size"=>$fsize, -cdata => $gene, style=>{'writing-mode'=>'tb','text-anchor'=>"middle",});		
		}
		$svg->rectangle(x=>$x0+($start-$startRef)/$unit,y=>$y0-$up,width=>($end-$start)/$unit, height=>$hei,style=>{'fill'=>$color,'stroke'=>$lgray, 'stroke-width'=>0,});
				
		$pos{$gene}[0] = [($x0 + ($start- $startRef) / $unit , $y0 - $up + $hei )];
        $pos{$gene}[1] = [($x0 + ($end  - $startRef) / $unit , $y0 - $up + $hei )];					
	}
	
	## draw other genomes
	foreach my $i (0..$#accs) {
	   next if ($i==0);
	   my $acc = $accs[$i];
	   my $y1 = $y0+$i*$top;
	   
	   	   
	   my $startAlt = $reg{$acc}[1];
	   my $endAlt = $reg{$acc}[2];
	   
	   $svg->text(id=>$acc, x=>$x0-30, y=>$y1, "font-family"=>"Arial","font-size"=>$fsize*2, -cdata => $acc, );	   
	   $svg->rectangle(x=>$x0,y=>$y1,width=>($endAlt-$startAlt)/$unit, height=>$h, style=>{'fill'=>$lgray,'stroke'=>$lgray, 'stroke-width'=>0,});
	}
	
	
	### draw the alignment
	#my @accs=("Col","An-1","C24","Cvi","Eri","Kyo","Ler","Sha");
	$x0 = 40;
	$y0 = 40;
	$top = 25;
	my $ya=$y0;
	foreach my $i (0..$#accs-1) {
	   last if ($flag!=2);
	   my $acc1 = $accs[$i];
	   my $acc2 = $accs[$i+1];
	   my @alns = @{$aln{$acc1}{$acc2}};
	   
	   my $start1 = $reg{$acc1}[1];
	   my $end1 = $reg{$acc1}[2];
	   
	   my $start2 = $reg{$acc2}[1];
	   my $end2 = $reg{$acc2}[2];
	   
	   my $ya = $top*$i + $y0; 
	   foreach my $blk (@alns) {
	   	  my @blk = @{$blk};
	   	  my ($st1,$en1,$st2,$en2) = @blk[1,2,4,5];
	   	  if ($blk[6] eq "-") {
	   	  	($st2,$en2) = ($en2, $st2);
	   	  }
	   	  $st1 = $start1 if ($st1 < $start1);
	   	  my $bx0 = $x0 +($st1 - $start1)/$unit; ##start point
	   	  my $by0 = $ya ;
	   	  $st2 = $start2 if ($st2 < $start2);
	   	  my $bx3 = $x0 + ($st2 - $start2)/$unit; ##end point 
	   	  my $by3 = $ya + $top  ;
	   	  
	   	  my ($bx1,$by1,$bx2,$by2); #two control points for Bezier curve
	   	  $bx1 = $bx0;
	   	  $by1 = ($by0+$by3)/2;
	   	  $bx2 = $bx3;
	   	  $by2 = $by1;
	   	 
	   	  my $curve = "M$bx0 $by0 C$bx1 $by1, $bx2 $by2, $bx3 $by3";
	   	  $en1 = $end1 if ($en1 > $end1);	   	  
	   	  $bx0 = $x0 + ($en1 - $start1)/$unit; ##start point
	   	  $by0 = $ya;
	   	  $en2 = $end2 if ($en2 > $end2);
	   	  $bx3 = $x0 + ($en2 - $start2)/$unit; ##end point 
	   	  $by3 = $ya + $top ;
	   	  	   	  
	   	  $bx1 = $bx0;
	   	  $by1 = ($by0+$by3)/2;
	   	  $bx2 = $bx3;
	   	  $by2 = $by1;
	   	  $curve = $curve . " L$bx3 $by3 C$bx2 $by2, $bx1 $by1, $bx0 $by0 z";
	   	  #next;
	   	  my $color=$vvlblue;
	   	  if ( $blk[-1] eq "SYN") {
	   	  	 $color = $lblue;
	   	  }elsif ($blk[-1] eq "INV") {
	   	  	$color=$lred;
	   	  }elsif ($blk[-1] eq "ITX" ) {
	   	  	$color=$lgreen;
	   	  }else {
	   	  	$color = $orange;
	   	  }
	   	  
	   	  my $tag = $svg->path(
             d  => $curve,             
		     style => { 'fill' => $color , "fill-opacity"=>"0.2"}
          );		   	     
	   }	   
    
	}
	
	
	## draw the ortholog
	foreach my $i (0..$#accs) {
	   next if ($i==0);
	   my $acc = $accs[$i];
	   my $y1 = $y0+$i*$top;
	   my $startAlt = $reg{$acc}[1];
	   my $endAlt = $reg{$acc}[2];
	   
	   foreach my $gene (sort keys %{$genes{$acc}}) {
	   	 my ($id,$start,$end,$direc) = @{$genes{$acc}{$gene}}; 	   	
	   	 my $color = $lblue;   
	   	 if ($rgenes{$acc}{$gene}) {
	   	 	if ($rgenes{$acc}{$gene}==1) {
			  	$color = $red; ##TNL,TN
			}elsif ($rgenes{$acc}{$gene}==2 ) {
			    $color = $orange; ##CNL,CN
			}elsif ($rgenes{$acc}{$gene}==3){
			    $color = $green; #other NBS-LRR
			}elsif ($rgenes{$acc}{$gene}==4){
			    $color = $gold; ##RLK
			}	elsif ($rgenes{$acc}{$gene}==5){
			    $color = $gray; ##RLP
			}else {
				$color = $puple; ##TM-CC
			}					
		 }
		 $svg->rectangle(x=>$x0+($start-$startAlt)/$unit,y=>$y1-$up,width=>($end-$start)/$unit, height=>$hei,style=>{'fill'=>$color,'stroke'=>$lgray, 'stroke-width'=>0,});
				 
		 ## draw the ortholog links	  
		 if ($ort{$gene}) {
		   my @orts = @{$ort{$gene}};
		   foreach my $ort (@orts) {
		   	 if (!$pos{$ort}) {
		   	 	print "No position: $gene\t$ort\n";
		   	 	next;
		   	 }
		   	 my @pos = @{$pos{$ort}};
		   	 my $bx0 = $pos[0][0]; ##start point
	   	     my $by0 = $pos[0][1];
	   	     
	   	     my $bx3 = $x0 + ($start - $startAlt) / $unit ; ##end point 
	   	     my $by3 = $y1 - $up ;
	   	     
	   	     my ($bx1,$by1,$bx2,$by2); #two control points for Bezier curve
	   	     $bx1 = $bx0;
	   	     $by1 = ($by0+$by3)/2;
	   	     $bx2 = $bx3;
	   	     $by2 = $by1;
	   	 
	   	     my $curve = "M$bx0 $by0 C$bx1 $by1, $bx2 $by2, $bx3 $by3";
		   	 
		   	 $bx0 = $pos[1][0]; ##start point
	   	     $by0 = $pos[1][1];
	   	     
	   	     $bx3 = $x0 + ($end   - $startAlt) / $unit ; ##end point 
	   	     $by3 = $y1 - $up ;
	   	  	   	  
	   	     $bx1 = $bx0;
	   	     $by1 = ($by0+$by3)/2;
	   	     $bx2 = $bx3;
	   	     $by2 = $by1;
	   	     $curve = $curve . " L$bx3 $by3 C$bx2 $by2, $bx1 $by1, $bx0 $by0 z";
		   	 
		   	 #$svg->line(x1=>($pos[0][0]+$pos[1][0])/2, y1=>($pos[0][1]+$pos[1][1])/2, 
		   	 #           x2=>$x0+($start-$startAlt)/$unit + ($end-$start)/(2*$unit), 
		   	 #           y2=>$y1 - $up, 
		   	 #           style=>{'fill'=>$red,'stroke'=>$gray,'stroke-width'=>0.5,});
		   	 #next;
		   	 next if (($rgenes{$acc}{$gene} or $rgenes{$accs[$i-1]}{$ort}) and ($flag!=1)) ; 
		   	 my $tag = $svg->path(
               d  => $curve,
               id => "$gene\-$ort",
		       style => { 'fill' => $llgray }
             );
		   	 
		     	
		   }
		   	
		 } 	    
	      $pos{$gene}[0] = [($x0 + ($start - $startAlt) / $unit  , $y1 - $up + $hei)]; # start
	      $pos{$gene}[1] = [($x0 + ($end   - $startAlt) / $unit  , $y1 - $up + $hei)]; 	# end 
	   }	      
	}
	
	
	
    open OUT,">$outSVG";
	print OUT $svg->xmlify;  
	#print "xxx $outSVG\n";
	
}



sub drawFig2 {
	# https://www.w3.org/TR/SVG/paths.html#PathDataLinetoCommands
	# https://metacpan.org/pod/Math::Bezier
	# https://metacpan.org/pod/SVG#path	
   my ($reg, $aln, $maxSize, $genes, $ort, $rgenes , $accs, $outSVG,$flag, $gr,$rgr) = @_;
   
   my @accs = @{$accs};
   my %aln = %{$aln};
   my %reg = %{$reg};
   	
   my %genes = %{$genes};
   my %ort = %{$ort};
   my %rgenes = %{$rgenes};
   my %gr = %{$gr};
   my %rgr = %{$rgr};
   
   my $n =0;
   my @colors = ("rgb(102,194,165)","rgb(252,141,98)","rgb(141,160,203)", "rgb(231,138,195)", "rgb(166,216,84)", "rgb(255,217,47)","rgb(229,196,148)","rgb(179,179,179)");
   foreach my $gr (sort keys %rgr) {
   	 if ($rgr{$gr}>0) {
   	 	$n+=1;
   	 	print "r gene group $gr\n";
   	 }
   }
   print "$outSVG : colored r gene group : $n\n";
   if ($n>8) {
   	 print "# >8 R gene groups \n";
   }
   
   
   
   my $unit=500; #scale	
   my $width = $maxSize/$unit + 120 ;
   my $height = ($#accs+1)*40 + 120;		
   my $svg= SVG->new(width=>$width,height=>$height);
   print "#figure size $outSVG $width  x  $height\n";     
	
	my $gold="rgb(255,215,0)";
	my $red="rgb(255,0,0)";
	my $lred="rgb(128,0,0)";
	my $blue="rgb(0,0,255)";
	my $lblue="rgb(0,0,128)";
	my $green="rgb(0,255,0)";
	my $lgreen="rgb(0,128,0)";
	my $gray="rgb(24,24,24)";
	my $llgray="rgb(196,196,196)";
	my $lgray="rgb(97,96,96)";
	my $black="rgb(0,0,0)";
	my $dgreen="rgb(53,128,0)";
    my $vvlblue = "rgb(239,243,255)";
    my $orange = "rgb(255, 127, 0 )";
    my $puple = "rgb(128,0,128)";
	##top fig for Col
	
	#draw the Col reg 
	my $x0 = 40;
	my $y0 = 40; 	
	my $top = 25; # distance between up and down regions
	my $hei = 10; #rectangle height for gene
	my $h = 2 ; ##  rectangle height for the region
	my $up = int(($hei-$h)/2); ## the height of up region of a gene
	my $fsize = 3 ; # font size  
	
    my $endRef = $reg{"Col"}[2];
	my $startRef = $reg{"Col"}[1];
	$svg->text(id=>"Col-0", x=>$x0-30, y=>$y0, "font-family"=>"Arial","font-size"=>$fsize*2, -cdata => "Col", );
	$svg->rectangle(x=>$x0,y=>$y0,width=>($endRef-$startRef)/$unit, height=>$h,style=>{'fill'=>$lgray,'stroke'=>$lgray, 'stroke-width'=>0,});
		
	#start position
	$svg->text(x=>$x0, y=>$y0-20, "font-family"=>"Arial","font-size"=>$fsize, -cdata => $startRef, style=>{'writing-mode'=>'tb','text-anchor'=>"middle",});
	#end position
	$svg->text(x=>$x0+($endRef-$startRef)/$unit, y=>$y0-20, "font-family"=>"Arial","font-size"=>$fsize, -cdata => $endRef, style=>{'writing-mode'=>'tb','text-anchor'=>"middle",});
	
	
    my %pos;
    
	foreach my $gene (sort keys %{$genes{"Col"}}) {
		my ($id,$start,$end,$direc) = @{$genes{"Col"}{$gene}};
		my $color = $lblue;      
		#print "gene $gene\n";  	
		if ($gr{$gene}) {
		   my $grNo = $gr{$gene};
		   if ($rgr{$grNo}>0) {
		  	$color = $colors[$rgr{$grNo}-1];
		  	print "$outSVG $grNo: $color \t $gene\n";
		  }				
							
			$svg->text(x=>$x0+($start-$startRef)/$unit, y=>$y0-20, "font-family"=>"Arial","font-size"=>$fsize, -cdata => $gene, style=>{'writing-mode'=>'tb','text-anchor'=>"middle",});		
		}else {
		  print "$gr{$gene}, $gene xxx\n";
		  if ($rgenes{"Col"}{$gene}) {
		  	$color=$gold;
		  }
		}
		$svg->rectangle(x=>$x0+($start-$startRef)/$unit,y=>$y0-$up,width=>($end-$start)/$unit, height=>$hei,style=>{'fill'=>$color,'stroke'=>$lgray, 'stroke-width'=>0,});
				
		$pos{$gene}[0] = [($x0 + ($start- $startRef) / $unit , $y0 - $up + $hei )];
        $pos{$gene}[1] = [($x0 + ($end  - $startRef) / $unit , $y0 - $up + $hei )];					
	}
	
	## draw other genomes
	foreach my $i (0..$#accs) {
	   next if ($i==0);
	   my $acc = $accs[$i];
	   my $y1 = $y0+$i*$top;
	   
	   	   
	   my $startAlt = $reg{$acc}[1];
	   my $endAlt = $reg{$acc}[2];
	   
	   $svg->text(id=>$acc, x=>$x0-30, y=>$y1, "font-family"=>"Arial","font-size"=>$fsize*2, -cdata => $acc, );	   
	   $svg->rectangle(x=>$x0,y=>$y1,width=>($endAlt-$startAlt)/$unit, height=>$h, style=>{'fill'=>$lgray,'stroke'=>$lgray, 'stroke-width'=>0,});
	}
	
	
	### draw the alignment
	#my @accs=("Col","An-1","C24","Cvi","Eri","Kyo","Ler","Sha");
	$x0 = 40;
	$y0 = 40;
	$top = 25;
	my $ya=$y0;
	foreach my $i (0..$#accs-1) {
	   last if ($flag!=2);
	   my $acc1 = $accs[$i];
	   my $acc2 = $accs[$i+1];
	   next if (!$aln{$acc1}{$acc2});
	   my @alns = @{$aln{$acc1}{$acc2}};
	   
	   my $start1 = $reg{$acc1}[1];
	   my $end1 = $reg{$acc1}[2];
	   
	   my $start2 = $reg{$acc2}[1];
	   my $end2 = $reg{$acc2}[2];
	   
	   my $ya = $top*$i + $y0; 
	   foreach my $blk (@alns) {
	   	  my @blk = @{$blk};
	   	  my ($st1,$en1,$st2,$en2) = @blk[1,2,4,5];
	   	  if ($blk[6] eq "-") {
	   	  	($st2,$en2) = ($en2, $st2);
	   	  }
	   	  $st1 = $start1 if ($st1 < $start1);
	   	  my $bx0 = $x0 +($st1 - $start1)/$unit; ##start point
	   	  my $by0 = $ya ;
	   	  $st2 = $start2 if ($st2 < $start2);
	   	  my $bx3 = $x0 + ($st2 - $start2)/$unit; ##end point 
	   	  my $by3 = $ya + $top  ;
	   	  
	   	  my ($bx1,$by1,$bx2,$by2); #two control points for Bezier curve
	   	  $bx1 = $bx0;
	   	  $by1 = ($by0+$by3)/2;
	   	  $bx2 = $bx3;
	   	  $by2 = $by1;
	   	 
	   	  my $curve = "M$bx0 $by0 C$bx1 $by1, $bx2 $by2, $bx3 $by3";
	   	  $en1 = $end1 if ($en1 > $end1);	   	  
	   	  $bx0 = $x0 + ($en1 - $start1)/$unit; ##start point
	   	  $by0 = $ya;
	   	  $en2 = $end2 if ($en2 > $end2);
	   	  $bx3 = $x0 + ($en2 - $start2)/$unit; ##end point 
	   	  $by3 = $ya + $top ;
	   	  	   	  
	   	  $bx1 = $bx0;
	   	  $by1 = ($by0+$by3)/2;
	   	  $bx2 = $bx3;
	   	  $by2 = $by1;
	   	  $curve = $curve . " L$bx3 $by3 C$bx2 $by2, $bx1 $by1, $bx0 $by0 z";
	   	  #next;
	   	  my $color=$vvlblue;
	   	  if ( $blk[-1] eq "SYN") {
	   	  	 $color = $lblue;
	   	  }elsif ($blk[-1] eq "INV") {
	   	  	$color=$lred;
	   	  }elsif ($blk[-1] eq "ITX" ) {
	   	  	$color=$lgreen;
	   	  }else {
	   	  	$color = $orange;
	   	  }
	   	  
	   	  my $tag = $svg->path(
             d  => $curve,             
		     style => { 'fill' => $color , "fill-opacity"=>"0.2"}
          );		   	     
	   }	   
    
	}
	
	
	## draw the ortholog
	foreach my $i (0..$#accs) {
	   next if ($i==0);
	   my $acc = $accs[$i];
	   my $y1 = $y0+$i*$top;
	   my $startAlt = $reg{$acc}[1];
	   my $endAlt = $reg{$acc}[2];
	   
	   foreach my $gene (sort keys %{$genes{$acc}}) {
	   	 my ($id,$start,$end,$direc) = @{$genes{$acc}{$gene}}; 	   	
	   	 my $color = $lblue;
	   	 if ($gr{$gene}) {
		   my $grNo = $gr{$gene};
		   if ($rgr{$grNo}>0) {
		  	 $color = $colors[$rgr{$grNo}-1];
		   }				
	   	 }
							
	   	 
		 $svg->rectangle(x=>$x0+($start-$startAlt)/$unit,y=>$y1-$up,width=>($end-$start)/$unit, height=>$hei,style=>{'fill'=>$color,'stroke'=>$lgray, 'stroke-width'=>0,});
				 
		 ## draw the ortholog links	  
		 if ($ort{$gene}) {
		   my @orts = @{$ort{$gene}};
		   foreach my $ort (@orts) {
		   	 if (!$pos{$ort}) {
		   	 	print "No position: $gene\t$ort\n";
		   	 	next;
		   	 }
		   	 my @pos = @{$pos{$ort}};
		   	 my $bx0 = $pos[0][0]; ##start point
	   	     my $by0 = $pos[0][1];
	   	     
	   	     my $bx3 = $x0 + ($start - $startAlt) / $unit ; ##end point 
	   	     my $by3 = $y1 - $up ;
	   	     
	   	     my ($bx1,$by1,$bx2,$by2); #two control points for Bezier curve
	   	     $bx1 = $bx0;
	   	     $by1 = ($by0+$by3)/2;
	   	     $bx2 = $bx3;
	   	     $by2 = $by1;
	   	 
	   	     my $curve = "M$bx0 $by0 C$bx1 $by1, $bx2 $by2, $bx3 $by3";
		   	 
		   	 $bx0 = $pos[1][0]; ##start point
	   	     $by0 = $pos[1][1];
	   	     
	   	     $bx3 = $x0 + ($end   - $startAlt) / $unit ; ##end point 
	   	     $by3 = $y1 - $up ;
	   	  	   	  
	   	     $bx1 = $bx0;
	   	     $by1 = ($by0+$by3)/2;
	   	     $bx2 = $bx3;
	   	     $by2 = $by1;
	   	     $curve = $curve . " L$bx3 $by3 C$bx2 $by2, $bx1 $by1, $bx0 $by0 z";
		   	 
		   	 #$svg->line(x1=>($pos[0][0]+$pos[1][0])/2, y1=>($pos[0][1]+$pos[1][1])/2, 
		   	 #           x2=>$x0+($start-$startAlt)/$unit + ($end-$start)/(2*$unit), 
		   	 #           y2=>$y1 - $up, 
		   	 #           style=>{'fill'=>$red,'stroke'=>$gray,'stroke-width'=>0.5,});
		   	 #next;
		   	 next if (($rgenes{$acc}{$gene} or $rgenes{$accs[$i-1]}{$ort}) and ($flag!=1)) ; 
		   	 next if ($gr{$gene} and $rgr{$gr{$gene}} > 0 );
		   	 my $tag = $svg->path(
               d  => $curve,
               id => "$gene\-$ort",
		       style => { 'fill' => $llgray }
             );
		   	 
		     	
		   }
		   	
		 } 	    
	      $pos{$gene}[0] = [($x0 + ($start - $startAlt) / $unit  , $y1 - $up + $hei)]; # start
	      $pos{$gene}[1] = [($x0 + ($end   - $startAlt) / $unit  , $y1 - $up + $hei)]; 	# end 
	   }	      
	}
	
	
	
    open OUT,">$outSVG";
	print OUT $svg->xmlify;  
	
	
}