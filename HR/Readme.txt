
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
