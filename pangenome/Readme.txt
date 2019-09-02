
pangenome can be built based on the whole genome sequence alignment or protein-coding genes ortholog clustering

# Pan-genome: genome sequence alignment

do all pairwise whole genome comparisons using MUMmer

prepare all chromosome length information in a file for each genome; e.g:
Chr1    30401407
Chr2    19417579
Chr3    23034411
Chr4    18785460
Chr5    26733864


python -u wga.pangenome.py -w ./pairwiseAssV2 -o ./ -g ../chrBed_v2/ &

# Pan-genome: protein-coding genes ortholog clustering
python pangenome.build.py -g AMPRIL.Alyrata.ortholog.groups.csv -o ./
