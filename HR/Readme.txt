
# synteny diversity
Synteny diversity was defined as the average fraction of non-syntenic sites found within all pairwise genome comparisons within
a given population. Here we denote synteny diversity as (Eq. 3) 

    𝜋𝑠𝑦𝑛=∑𝑥𝑖𝑥𝑗𝜋𝑖𝑗, (3) 

where xi and xj refer to the frequencies of sequence i and j and πij to the average probability of a position to be non-syntenic
between sequence i and j . Note, πsyn can be calculated in a given region or for the entire genome. However even when calculated
for small regions the annotation of synteny still needs to be established within the context of the whole genomes to avoid false
assignments of homologous but non-allelic sequence. Here we used the annotation of SyRI to define syntenic regions. πsyn values
can range from 0 to 1, with higher values referring to a higher average degree of non-syntenic regions between the genomes.


# step 1: pairwise whole genome alignment and synteny identification 
Before caculate synteny diversity, please run all pairwise whole genome comparison using the tool nucmer in the MUMmer4 package
and run SyRi to identify the syntenic and rearranged regions for each comparison.

Let's assume all the alignments of between our 8 A.thaliana genomes in a folder like below
/xxx/pairwiseWGA
/xxx/pairwiseWGA/An-1
/xxx/pairwiseWGA/An-1/C24
/xxx/pairwiseWGA/An-1/Cvi
...
/xxx/pairwiseWGA/An-1/Sha
...
...
...
/xxx/pairwiseWGA/Sha
/xxx/pairwiseWGA/Sha/An-1
/xxx/pairwiseWGA/Sha/C24
...

# step 2:  get the coordinates of syntenic regions in all genomes
  bedtools multiinter -i /xxx/pairwiseWGA/Col/*/*.wga.syn.block.txt  -names An-1 C24 Cvi Eri Kyo Ler Sha >Col.syn.txt
  
  Note: file format for seven *.wga.syn.block.txt files at below. Tab seprated ,
 
  Chr1    1926    23160   Chr1    2148    23388   +   AB  SYN
  Chr1    23138   44997   Chr1    23454   45322   +   AB  SYN
  Chr1    45229   68996   Chr1    45320   69077   +   AB  SYN
  Chr1    69082   76826   Chr1    70212   77969   +   AB  SYN
  Chr1    77557   132220  Chr1    77966   132688  +   AB  SYN

  The latest version of SyRI does not generate outputs like the above, while it can be easily generated from the SyRI output. e.g: grep SYNAL syri.out|cut -f 1-3,6-8,10  The 7th and 8th columns are not necessary. 
  
  for k in {An-1,C24,Cvi,Eri,Kyo,Ler,Sha}; do
    for j in {1..5}; do show-aligns -r /xxx/pairwiseWGA/Col/$k/out_m_i90_l100.delta Chr$j Chr$j >/xxx/pairwiseWGA/Col/$k/out_m_i90_l100.chr$j.aligns ; done ;
  done
  
  for k in {An-1,C24,Cvi,Eri,Kyo,Ler,Sha};  do cat $k/out_m_i90_l100.chr*.aligns >$k/$k.aligns ;done &
  
  awk '{if ($4==7) print}' Col.syn.txt >Col.syn.all.txt  ##change the "7" according to the number of your genomes
  
  perl get.all.syn.coord.pl ./Col.syn.all.txt ../../results/pairwiseAssV2/Col/ ./Col.syn.all.coords.txt &  ##maybe, please modify the line 86 and 98 accordingly if you have different number of genomes
  
  The output file Col.syn.all.coords.txt :
  $head -n3 syn/Col.syn.all.coords.txt (every three columns reprent a region of one genome)
    Chr1	1926	17011	An-1	2148	17238	C24	849	15976	Cvi	791	15861	Eri	2276	17379	Kyo	813	15903	Ler	842	15933	Sha	2946	18118
    Chr1	18731	44997	An-1	18941	45322	C24	17696	43984	Cvi	17564	43835	Eri	17379	43675	Kyo	17623	43884	Ler	17653	43904	Sha	18118	44406
    Chr1	45229	55676	An-1	45320	55772	C24	44216	54685	Cvi	44068	54515	Eri	43907	54352	Kyo	44116	54563	Ler	44136	54583	Sha	44638	55092
    
   $head Col.syn.all.coords.txt2 (add three columns for Col-0 <reference> from the above result file Col.syn.all.txt, compared to Col.syn.all.coords.txt)
Chr1	1926	17011	An-1	2148	17238	C24	849	15976	Col	1926	17011	Cvi	791	15861	Eri	2276	17379	Kyo	813	15903	Ler	842	15933	Sha	2946	18118
Chr1	18731	44997	An-1	18941	45322	C24	17696	43984	Col	18731	44997	Cvi	17564	43835	Eri	17379	43675	Kyo	17623	43884	Ler	17653	43904	Sha	18118	44406
Chr1	45229	55676	An-1	45320	55772	C24	44216	54685	Col	45229	55676	Cvi	44068	54515	Eri	43907	54352	Kyo	44116	54563	Ler	44136	54583	Sha	44638	55092

    
# step 3: caculate synteny diversity for every postion of the genome
  perl calculate.syn.diversity.pl ./Col.syn.all.coords.txt2 /xxx/pairwiseAssV2/ chrBed_v2 ./syn.diversity.position.Col.txt 

Note: the folder chrBed_v2 contains xx.leng.txt (format: chromosome\tlength\n) files for each genome e.g: $cat chrBed_v2/An-1.leng.txt
Chr1	30401407
Chr2	19417579
Chr3	23034411
Chr4	18785460
Chr5	26733864

Alternatively, in order to avoid too much memory occupation, caculate each pairwise per chromosome, then merge and calculate
nohup perl ../../scripts/run.cal.syn.div.pl ./Col.syn.all.coords.txt2 ../../results/pairwiseAssV2/ ../../chrBed_v2/ ../../scripts/calculate.syn.diversity.pairwise.chr.pl ./splitChr/ > np.chr4 &

#
nohup perl ../../scripts/syn.div.merge.pl ./splitChr/ Chr2 ../../chrBed_v2/Col.leng.txt ./splitChr/Chr2.syn.div.pos.txt > np.merge2


# step 4: caculate synteny diversity in a sliding window
  for k in {1..5}; do   perl calculate.syn.diversity.window.pl ./splitChr/Chr$k.syn.div.pos.txt 5000 1000 splitChr/Chr$k.syn.div.win5kb.step1kb.txt & done &
  cat syn.div.win5kb.step1kb.txt  > 
  
  Note: 
    the result file syn.div.win5kb.step1kb.txt, the last column shows the average syntenty diversity in a window
      Chr1	1	5000	0	0.220
      Chr1	1001	6000	1	0.090
      Chr1	2001	7000	2	0.000
      Chr1	3001	8000	3	0.000
      Chr1	4001	9000	4	0.000
      

# other downstream analysis
## find Hotspot of Rearrangements (HOR)
  awk '{if ($5>0.5)print}' syn.div.win5kb.step1kb.txt |bedtools merge -i - -d 2000 >syn.div.win5kb.step1kb.HDR.bed
  Note: 
      1. here we set the cutoff value 0.5, suggesting more than two divergent "synteny haplotypes" existing in our eight genomes.
      2. we also merge ajacent HORs with distance shorter than 2kb.
           
##  visualize the gene arrangments in all HOR(HDR)
 perl HDR.gene.scheme.2.pl Col.syn.all.coords.txt2 ./HDR.clu.bed /xxx/pairwiseAssV2/ AMPRIL.ortholog.groups.csv geneBed2/ Rgenes 50000 ./cluWin50kb2
 
 Note: 
  1. $head -n3 file HDR.clu.bed
  Chr1	1441001	1447000	A	A	Chr1_1441001_1447000
  Chr1	2293001	2298000	A	A	Chr1_2293001_2298000
  Chr1	2788001	2801000	A	A	Chr1_2788001_2801000
  2. gene family clustering result file "AMPRIL.ortholog.groups.csv" , generated by the tool orthofinder, e.g:
  OG0001305:	ATAN-1G83000,ATAN-4G46710	ATC24-1G78390,ATC24-4G49420	AT1G65480,AT4G20370	ATCVI-1G79370,ATCVI-4G43930	ATERI-1G76900,ATERI-4G45410	ATKYO-1G78110,ATKYO-4G46310	ATLER-1G78820,ATLER-4G47010	ATSHA-1G79630
  3. the folder geneBed2 contains 8 xxx.protein-coding.genes.v2.5.bed files for the coordinate of protein-coding genes 
      Chr1	4779	6663	+	ID=ATSHA-1G10010;Note=protein_coding_gene
      Chr1	7937	9721	-	ID=ATSHA-1G10020;Note=protein_coding_gene
      Chr1	12952	14028	-	ID=ATSHA-1G10030;Note=protein_coding_gene

  4. the folder Rgenes contains R gene annotation result file for all eight genomes . e.g:
    chr1	194627	197030	+	ATSHA-1G10570	RLK	RLK
    chr1	890057	891985	-	ATSHA-1G12860	TM-CC	TM-CC
    chr1	1721562	1725011	+	ATSHA-1G15300	RLK	RLK

##  visualize the gene arrangments all R gene clusters
    perl R.gene.cluster.wga.ortho.scheme.pl Col.syn.all.coords.txt2 ./R.gene.cluster.bed AMPRIL.ortholog.groups.csv geneBed2/ Rgenes 50000 ./cluWin50kb 
