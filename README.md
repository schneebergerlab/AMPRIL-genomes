# Overview
These scripts are used for the project of seven Arabidopsis thaliana genomes assembly, annotation and comparisons.


# System Requirements
The scripts have been tested on Linux (4.4.0-97-generic #120-Ubuntu)
Perl (5.22.2)
Python (2.7.15)
R (3.5.2)

# Installation
git clone https://github.com/schneebergerlab/AMPRIL-genomes.git
All scripts can be run directly.

# Usage

## workflow for protein-coding gene annotation 
### step 1: run pipeline for protein-coding gene annotation
  
  python ../scripts/evm.pasa.integrate.pipeline.py -f ./annotation.config
  This script will do 
  
  1) protein sequence alignment using exonerate  
  2) RNAseq reads mapping using hisat2		
  3) ab initio prediction using AUGUSTUS, GlimmerHMM, SNAP		
  4) merge results from 1),2) and 3)
  5) run EVM to integrate the results to get consensus gene models
  6) get the gene, protein, CDS sequences based on the gene model gff3 file and reference fasta file
  The results will be included in the file evm.all.gff3
  
### step 2: repeat annotation
  RepeatMasker -species arabidopsis -gff -dir ./ -pa 20 ../../reference/chr.all.v2.0.fasta
  perl ../../../scripts/repeat.classfied.gff3.pl ./chr.all.v2.0.fasta.out.gff ./chr.all.v2.0.fasta.out ./chr.all.v2.0.fasta.repeats.ann.gff3 repeat.ann.stats &
  egrep -v 'Low|Simple|RNA|other|Satellite' chr.all.v2.0.fasta.repeats.ann.gff3 |cut -f 1,4,5,9 >chr.all.v2.0.TE.bed

### step 3: identify TE-related genes
  perl ../../../scripts/remove.TErelated.genes.pl ../../EVM_PASA/evm.annotation.protein.fasta ../../EVM_PASA/evm.annotation.gene.fasta ../RepeatMasker/chr.all.v2.0.TE.bed ../../EVM_PASA/evm.all.gff3 ./ 

### step 4: non-coding gene annotation
  cmscan --cpu 20 --tblout Rfam.scan.out ../../../data/Rfam/Rfam.cm ../chr.all.v2.0.split.fasta
  perl  ../../../scripts/noncoding.infernal.output.parser.pl ./Rfam.scan.out ./ >nohup.out & 

### step 5: get the initial version
   nohup python ./scripts/gene.id.update.py -i ./ -v 2.0 >log/gene.id.update.log &
   
### step 6: evaluation
  awk '{if ($3!="miRNA_primary_transcript" && $3!="pseudogenic_exon" && $3!="pseudogenic_transcript" && $3!="pseudogenic_tRNA" && $3!="transposon_fragment" && $3!="mRNA" && $3!="protein" && $3 !="CDS" && $3!="exon" && $3!="five_prime_UTR" && $3!="three_prime_UTR" && $3!="lnc_RNA") print}' Araport11_GFF3_genes_transposons.201606.gff |grep -v '^#' |cut -f 1,4,5,9 |grep -vP '\.\d+\;Parent' |grep -v  'ChrC|ChrM' >Araport11_gene.TE.chr1-5.bed &

cd /AMPRILdenovo/annotation/Sha/evaluation
mkdir blastnCol misannotation update 

#### 1) gene seq. blastn
cd blastnCol
ln -s ../../version/Sha.v2.0.gene.fasta ./gene.fasta
blastn -query ./gene.fasta  -db TAIR10_chr_all.fas -num_threads 20 -evalue 1e-5 -out gene.blastn.Col.out &

perl annotation.evaluation.by.geneBlastnCol.pl ./gene.blastn.Col.out
awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' ../blastnCol/gene.blastn.besthit.out |sed 's/chloroplast/ChrC/g' |sed 's/mitochondria/ChrM/g' |awk '{if (!/Chr/) print "Chr"$0;else print}' |sed 's/ChrNone/None/g' >query.gene.blastn.besthit.out

blastn -query Araport11.prot.genomic.seq.fasta -db chr.all.v2.0.fasta -out arap11.prot-gene.blastn.out -evalue 1e-5 &
perl ../../../scripts/assembly.eval.arap11.single.pl ./arap11.prot-gene.blastn.out ../../../../tair10/Araport11/Araport11.prot.genomic.seq.fasta ./araport11.gene.blastn.assembly.out &

#### 2) prot seq. blastp and gene family clustering
cd /AMPRILdenovo/genefamily/blastpAraport11/Sha
mkdir prot orthomcl; awk '{if (/>/) print $1;else print}' ../../../annotation/Sha/version/Sha.1.0.protein.fasta  |sed 's/>/>Sha|/g' > prot/Sha.fasta
cd prot/; ln -s ../../Col.fasta ./ ; cd ../ ; cat prot/*.fasta > Col-Sha.fasta ; makeblastdb -in ./Col-Sha.fasta -dbtype prot ;sed 's/Kyo/Sha/g' ../Kyo/orthomcl/orthomcl.config  >orthomcl/orthomcl.config; 

blastp -query Col-Sha.fasta -db Col-Sha.fasta -num_threads 40 -evalue 1e-10 -outfmt 6 -out blastout ; orthomclInstallSchema ./orthomcl/orthomcl.config install.sql.log; grep -P "^[^#]" blastout > blastresult; orthomclBlastParser blastresult prot > orthomcl/similarSequences.txt; perl -p -i -e 's/\t(\w+)(\|.*)orthomcl/\t$1$2$1/' orthomcl/similarSequences.txt; perl -p -i -e 's/0\t0/1\t-181/' orthomcl/similarSequences.txt; cd orthomcl ; orthomclLoadBlast ./orthomcl.config similarSequences.txt ; orthomclPairs ./orthomcl.config orthomcl_pairs.log cleanup=all  ; orthomclPairs ./orthomcl.config orthomcl_pairs.log cleanup=no; orthomclDumpPairsFiles ./orthomcl.config ;  mcl mclInput --abc -I 1.5 -o mclOutput -te 20; orthomclMclToGroups group 1 < mclOutput > groups.txt &


#### 3) find mis-merging, mis-spliting, missing, mis-annotated genes
/AMPRILdenovo/annotation/Cvi/evaluation
cd misannotation

ln -s ../../../../tair10/Araport11/Araport11_gene.bed Araport11.protein.bed
ln -s ../../../../genefamily/blastpAraport11/Col.fasta Araport11.protein.fasta

cd ../../version/ ; awk '{if ($3=="gene") print}' C24.protein-coding.genes.v1.0.gff |cut -f 1,4,5,7,9 >C24.protein-coding.genes.v1.0.bed ; 
ln -s ../../version/Sha.protein-coding.genes.v1.0.bed query.protein.bed
ln -s ../../../../genefamily/blastpAraport11/Sha/prot/Sha.fasta query.protein.fasta

ln -s  ../../../../genefamily/blastpAraport11/Sha/blastresult
ln -s  ../../../../genefamily/blastpAraport11/Sha/orthomcl/groups.txt ./

awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' ../../../../assembly/ShaNew/evaluation/Araport11blastn/araport11.gene.besthit.out >Araport11.gene.blastn.besthit.out
awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' ../blastnCol/gene.blastn.besthit.out >query.gene.blastn.besthit.out

##input files: 
 >>> Araport11 and accession protein region bed files and sequences files.
   grep gene  ../../repeat/TErelated/annotation.genes.gff|cut -f 1,4,5,9 |sed 's/TU/model/g' >query.prot.gene.bed   
 >>> Blastp result of accession proteins against Araport11 proteins
 >>> OrthoMCL clustering result between accession and Araport11 proteins
     grep AT  ../../../../genefamily/AssV2tmp/An-1/prot/Results_Aug11/Orthogroups.txt |grep evm >groups.txt
 >>> Blastn result of Araport11 gene sequences against the accession assembly (Blastn result 1)
     awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' Col.prot.besthit.out |grep -P 'AT\d' >Col.prot.besthit.out2
 >>> Blastn result of accession gene sequences against Col-0 genome sequences. (Blastn result 2)
     awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' query.prot.besthit.out >query.prot.besthit.out2

##output files:
 >>> potential.mis-merged.gene.txt (function: findMisMer, based on the result of blastp between Col and the Accession)
 >>> potential.mis-spliting.gene.txt (function: findMisSplit. based on the result of blastp between Col and the Accession)	
 >>> potential.query.un-assembled.gene.txt (findMissingMisann. based on the result of Col gene blastn against the accession's assembly)
 >>> potential.missing.gene.txt (findMissingMisann. based on the result of Col gene blastn against the accession's assembly)
 >>> potential.mis-exon-intron.gene.txt
 >>> potential.mis-split.gene.by.blastn.txt
 >>> potential.mis-merge.gene.by.blastn.txt
 >>> potential.m-vs-m.toBeChecked.by.blastn.txt
 >>> futher.check.list
 >>> potential.mis-annotated.gene.txt
 >>> Araport11.ungrouped.gene.analysis.stat
 >>> Araport11.ungrouped.gene.analysis.txt
 >>> Acc(xxx).ungrouped.gene.analysis.stat
 >>> Acc(xxx).ungrouped.gene.analysis.txt
 >>> xxx.genes.to.be.updated.added.txt (based on the potential.xxx.txt from blastn-based analysis)
 >>> xxx.genes.to.be.updated.added.srt.txt	

##to find genes:
 >>> mis-merged (Blastp result, blastn)
 >>> mis-split (Blastp result, blastn)
 >>> wrong exon-intron structure (blastn)
 >>> false protein-coding genes (not annotated in Araport 11 but actually they were assembled in the Col-0 genome) (Blastn result 2)
 >>> missing genes (not annotated in accession, but actually they were assembled)  (Blastn result 1)

   
nohup python -u ../../../scripts/annotation.evaluate.find-mis.py -g ./groups.txt -o ./run2 -n  Col.prot.besthit.out2 -c query.prot.besthit.out2 -p blastp.result -s Col.prot.gene.bed -q query.prot.gene.bed -x Col.prot.fasta -y query.prot.fasta -a Col.gene.LoF.txt -b query.gene.LoF.txt -r ../../RNAseq/hisat2/rnaseq.4evm.gff >np.run2&

### step 7: update    
#### 1) run scipio to align the protein sequences from Araport11 annoation
  python ../../../scripts/run.scipio.py -i ../../../data/protein/Araport11-split -o ./  -r ./reference/chr.all.v2.0.fasta >run.log &
  
  ln -s ../misannotation/Cvi.genes.to.be.updated.added.srt.txt ./; grep chr Cvi.genes.to.be.updated.added.srt.txt >genes.to.be.updated.txt;  grep -v chr Cvi.genes.to.be.updated.added.srt.txt >>genes.to.be.updated.txt ; ln -s ../../abinitio/abinitio.4evm.gff ; cp ../../../An-1/evaluation/update2/srt.gff.pl ./; perl ./srt.gff.pl ./abinitio.4evm.gff SNAP SNAP.ann.gff; perl srt.gff.pl ./abinitio.4evm.gff GlimmerHMM GlimmerHMM.ann.gff  ; cat ../../augustus/genome.chunk.*.gff |egrep -v '#|intron|transcription|codon' > ../../augustus/augustus.ann.gff; ln -s ../../augustus/augustus.ann.gff; ln -s ../misannotation/Araport11.gene.blastn.besthit.out2 ./ ; ln -s ../../../../tair10/Araport11/Araport11_genes.201606.pep.repr.fasta; ln -s ../misannotation/Araport11.protein.bed ./;  ln -s ../../version/Cvi.protein-coding.genes.v1.0.gff ./ ; cat ../../scipio/run2/splitOut/protein.chunk.*.gff >../../scipio/run2/splitOut/scipio.gff; ln -s  ../../scipio/run2/splitOut/scipio.gff ./ ;  ln -s ../../../../wga/results/Cvi/Cvi.wga.snp.indel.gene.LOF.txt ./wga.snp.indel.gene.LoF.txt

#### 2) update
	input files:
		gene.to.be.updated.txt
		xx.protein-coding.genes.v1.0.gff
		scipio.gff
		wga.snp.indel.gene.LoF.txt
		augustus.ann.gff; SNPA.ann.gff; GlimmerHMM.ann.gff
		Col gene blastn best hit out
		Col Araport11 protein sequence and bed files
		ChrCM.txt (a few of organella contigs were not removed in the assembly process)
	output files:
		genes.to.be.updated.txt2
		updated.gff
		updated.rmdup.gff
		updated.highConf.gff
		updated.highConf.prot.fasta
	Method:
		1) check the LoF information resulting from WGA-based SNPs and InDel annotation, and the Col protein sequence alignment result from Scipio, add the update information :
			ChrCM: ChrC or ChrM genes
			LowConf: low confident genes
			unchange: keep the previous annotation, 
			ChangeSci: annotate based on Scipio result (checkScipio: start codon, stop codon, splice sites, frame-shift, premature stop-codon gain; check AugGenes snapGenes glimGene, check AugGenes2[ab initio], check GeneWise)
			ChangeAug: annotate based on Augustus ( check AugGenes snapGenes glimGene, check AugGenes2[ab initio], check GeneWise)
			ChangeSciAug: annotate based on scipio and augustus (checkScipio, check AugGenes snapGenes glimGene, check AugGenes2[ab initio], check GeneWise)
			not-add: not annotated	
	
		2) prepare the other annotation result from Augustus-evidence-based method, SNAP ab initio and GlimmerHMM ab initio result
	 		cat ../../abinitio/augustus/augustus.hint.chunk.*.gff |egrep -v '#|intron|transcription|codon' > ../../augustus/augustus.ann.gff
		3) update 
  
   nohup python -u ../../../scripts/update.misann.genes.py -u genes.to.be.updated.txt -g annotation.genes.gff -o ./run2 -s scipio.gff -x Col.gene.LoF.txt -y query.gene.LoF.txt -c ChrCM.txt -a augustus.ann.gff -n SNAP.4evm.gff -l glimmerhmm.4evm.gff -b ./Col.gene.blastn.besthit.bed  -f ../../reference/chr.all.v2.0.fasta -p Col.prot.fasta -i ./Col.prot.gene.bed >update.run2.log  &
  
  python ../../scripts/annotation.gene.ID.update.py -i update2/updated.highConf.gff -n ../version/An-1.genes.annotation.v2.0.gff -o ../version -v v2.5 -a An-1 -g ../reference/chr.all.v2.0.fasta &

 
# Workflow for pan-genome analysis

pangenome can be built based on the whole genome sequence alignment or protein-coding genes ortholog clustering

## Pan-genome: genome sequence alignment

do all pairwise whole genome comparisons using MUMmer

prepare all chromosome length information in a file for each genome; e.g:
Chr1    30401407
Chr2    19417579
Chr3    23034411
Chr4    18785460
Chr5    26733864

python -u wga.pangenome.py -w ./pairwiseAssV2 -o ./ -g ../chrBed_v2/ &

## Pan-genome: protein-coding genes ortholog clustering
python pangenome.build.py -g AMPRIL.Alyrata.ortholog.groups.csv -o ./


# Workflow for Structural variations calling for each assembled genome.
All assemblies were aligned to the reference sequence (TAIR10) using nucmer from the MUMmer4 toolbox with parameter setting “-max -l 40 -g 90 -b 100 -c 200”. The resulting alignments were further filtered for alignment length (>100) and identity (>90). Structural rearrangements and local variations were identified using SyRI (https://github.com/schneebergerlab/syri).

For details of the file format, please check https://schneebergerlab.github.io/syri/fileformat.html

Example of commands used:
  nucmer --maxmatch  -l 40 -g 90 -c 100 -b 200 -t 20 Col.fasta An-1.fasta
  delta-filter -m -i 90 -l 100 out.delta > out_m_i90_l100.delta
  show-coords -THrd out_m_i90_l100.delta > out_m_i90_l100.coords
  syri -c out.chrom.coords -d out_m_i90_l100.delta -r Col.fasta -q An-1.fasta --nc 5 --all -k


# Workflow for the analysis of synteny diversity

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

## get coordinates of syntenic regions in all genomes
  bedtools multiinter -i ../../results/pairwiseAssV2/Col/*/*.wga.syn.block.txt  -names An-1 C24 Cvi Eri Kyo Ler Sha >Col.syn.txt
  for k in {1..5};do show-aligns -r out_m_i90_l100.delta Chr$k chr$k >out_m_i90_l100.chr$k.aligns ;done &
  for k in {An-1,C24,Cvi,Eri,Kyo,Ler,Sha};do cat $k/out_m_i90_l100.chr*.aligns >$k/$k.aligns ;done &
  bedtools multiinter -i ../../results/pairwiseAssV2/Col/*/*.wga.syn.block.txt  -names An-1 C24 Cvi Eri Kyo Ler Sha >Col.syn.txt
  perl ../../scripts/get.all.syn.coord.pl ./Col.syn.all.txt ../../results/pairwiseAssV2/Col/ ./Col.syn.all.coords.txt &

## caculate synteny diversity for every postion of the genome
perl ../../scripts/calculate.syn.diversity.pl ./Col.syn.all.coords.txt2 ../../results/pairwiseAssV2/ ../../chrBed_v2 ./syn.diversity.position.Col.txt 

## caculate synteny diversity in a sliding window
for k in {1..5}; do   perl ../../scripts/calculate.syn.diversity.window.pl ./splitChr/Chr$k.syn.div.pos.txt 5000 1000 splitChr/Chr$k.syn.div.win50kb.step5kb.txt & done &

### find HOR (HDR)
awk '{if ($5>0.5)print}' syn.div.win5kb.step1kb.txt |bedtools merge -i - -d 2002 |bedtools intersect -a - -b ../../../tair10/centromere_Giraut2011.bed -wao |awk '{if ($7==0) print $0"\tA";else print $0"\tC"}' |cut -f 1-3,8  >syn.div.win5kb.step1kb.HDR.bed

### gene arrangment in the HOR(HDR)
 nohup perl ../../scripts/syndiv/HDR.gene.scheme.2.pl ../00_synDiv/Col.syn.all.coords.txt2 ./HDR.clu.bed ../../01_syri/pairwiseAssV2/ ../../../genefamily/AMPRIL/ver3/AMPRIL.ortholog.groups.csv ../../../genefamily/AMPRIL/ver3/geneBed2/ ../../../genefamily/AMPRIL/ver3/Rgenes/ann/ 50000 ./cluWin50kb2 > np.log2 &

### gene arrangment in the HOR (HDR)
nohup perl ../../../scripts/Rgenes/R.gene.cluster.wga.ortho.scheme.pl ../../../../wga/07_synDiversity/00_synDiv/Col.syn.all.coords.txt2 ./R.gene.cluster.bed ../../../../wga/01_syri/pairwiseAssV2/ ../AMPRIL.ortholog.groups.csv ../geneBed2/ ./ann/ 20000 ./tmp >tmp.log&

  

# Citation:
Chromosome-level assemblies of multiple Arabidopsis genomes reveal hotspots of rearrangements with altered evolutionary dynamics
Wen-Biao Jiao, Korbinian Schneeberger. bioRxiv 738880; doi: https://doi.org/10.1101/738880

