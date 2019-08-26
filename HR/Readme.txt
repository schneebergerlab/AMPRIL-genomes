
synteny diversity

Before caculate synteny diversity, please run all pairwise whole genome comparison using MUMmer and run SyRi to identify the syntenic and rearranged regions for each comparison.

Let's assume all the alignments in a folder like below

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

# get coordinates of syntenic regions in all genomes
  bedtools multiinter -i ../../results/pairwiseAssV2/Col/*/*.wga.syn.block.txt  -names An-1 C24 Cvi Eri Kyo Ler Sha >Col.syn.txt
  for k in {1..5};do show-aligns -r out_m_i90_l100.delta Chr$k chr$k >out_m_i90_l100.chr$k.aligns ;done &
  for k in {An-1,C24,Cvi,Eri,Kyo,Ler,Sha};do cat $k/out_m_i90_l100.chr*.aligns >$k/$k.aligns ;done &
  bedtools multiinter -i ../../results/pairwiseAssV2/Col/*/*.wga.syn.block.txt  -names An-1 C24 Cvi Eri Kyo Ler Sha >Col.syn.txt
  perl ../../scripts/get.all.syn.coord.pl ./Col.syn.all.txt ../../results/pairwiseAssV2/Col/ ./Col.syn.all.coords.txt &

# caculate synteny diversity for every postion of the genome
perl ../../scripts/calculate.syn.diversity.pl ./Col.syn.all.coords.txt2 ../../results/pairwiseAssV2/ ../../chrBed_v2 ./syn.diversity.position.Col.txt 

# caculate synteny diversity in a sliding window
for k in {1..5}; do   perl ../../scripts/calculate.syn.diversity.window.pl ./splitChr/Chr$k.syn.div.pos.txt 5000 1000 splitChr/Chr$k.syn.div.win50kb.step5kb.txt & done &

## find HR (HDR)
awk '{if ($5>0.5)print}' syn.div.win5kb.step1kb.txt |bedtools merge -i - -d 2002 |bedtools intersect -a - -b ../../../tair10/centromere_Giraut2011.bed -wao |awk '{if ($7==0) print $0"\tA";else print $0"\tC"}' |cut -f 1-3,8  >syn.div.win5kb.step1kb.HDR.bed

## gene arrangment in the HR(HDR)
 nohup perl ../../scripts/syndiv/HDR.gene.scheme.2.pl ../00_synDiv/Col.syn.all.coords.txt2 ./HDR.clu.bed ../../01_syri/pairwiseAssV2/ ../../../genefamily/AMPRIL/ver3/AMPRIL.ortholog.groups.csv ../../../genefamily/AMPRIL/ver3/geneBed2/ ../../../genefamily/AMPRIL/ver3/Rgenes/ann/ 50000 ./cluWin50kb2 > np.log2 &

##gene arrangment in the HR(HDR)
nohup perl ../../../scripts/Rgenes/R.gene.cluster.wga.ortho.scheme.pl ../../../../wga/07_synDiversity/00_synDiv/Col.syn.all.coords.txt2 ./R.gene.cluster.bed ../../../../wga/01_syri/pairwiseAssV2/ ../AMPRIL.ortholog.groups.csv ../geneBed2/ ./ann/ 20000 ./tmp >tmp.log&
