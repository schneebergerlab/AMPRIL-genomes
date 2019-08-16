
pangenome can be built based on the whole genome sequence alignment or protein-coding genes ortholog clustering
#Method 1: genome sequence alignment
do all pairwise whole genome comparisons using MUMmer

python -u wga.pangenome.py -w ./pairwiseAssV2 -o ./ -g ../chrBed_v2/ &
#Method 2: protein-coding genes ortholog clustering
python pangenome.build.py -g AMPRIL.Alyrata.ortholog.groups.csv -o ./
