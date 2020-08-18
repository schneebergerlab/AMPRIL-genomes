
pangenome can be built based on the whole genome sequence alignment or protein-coding genes ortholog clustering

# Pan-genome based on genome sequence alignment

## step 1: do all pairwise whole genome comparisons using the tool nucmer in the package MUMmer

Let's assume all the alignments of between our 8 A.thaliana genomes in a folder like below
/xxx/pairwiseWGA
/xxx/pairwiseWGA/An-1
/xxx/pairwiseWGA/An-1/C24 (files: out_m_i90_l100.coords, C24.wga.block.txt, C24.del.bed, C24.ins.bed)
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

Note:
1. the alignment ".coords" result file, e.g out_m_i90_l100.coords:
  1	2251	14107439	14105198	2251	2242	93.74	1	-1	000041F	chr3
  2253	5758	11159799	11163301	3506	3503	90.71	1	1	000041F	chr5

>  prepare the bed formated file "xx.wga.block.txt" from the .coords file, e.g: C24.wga.block.txt
  Chr1	1	503	Chr4	19336140	19335640
  Chr1	751	1354	Chr1	626	2	
  Chr1	2148	4526	Chr1	849	3227	
  Chr1	4513	23388	Chr1	3268	22113	
> prepare the bed formated file "xx.del.bed" and "xx.ins.bed"
  bedtools complement -i -g C24.
  
  
## step 2: get pan-genome
python -u wga.pangenome.py -w ./pairwiseAssV2 -o ./ -g ../chrBed_v2/ &

Note: the folder chrBed_v2 contains xx.leng.txt (format: chromosome\tlength\n) files for each genome 
e.g: $cat chrBed_v2/An-1.leng.txt
Chr1	30401407
Chr2	19417579
Chr3	23034411
Chr4	18785460
Chr5	26733864



# Pan-genome: protein-coding genes ortholog clustering
python pangenome.build.py -g AMPRIL.ortholog.groups.csv -o ./
