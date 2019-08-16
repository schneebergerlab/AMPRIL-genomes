
bsub206 blastn -query ../../../../tair10/TAIR10_seq_20110103_representative_gene_model_updated -db ../../version/chr.all.v1.0.fasta -out tair10.gene.blastn.assembly.out -evalue 1e-5 -num_threads 6 &

## Reference :
#(1) Cheng_et_al-2017-The_Plant_Journal Araport11: a complete reannotation of the Arabidopsis thaliana reference genome

##1 protein-coding gene annotation
#protein : TAIR10/Araport11

## RNAseq preprocessing
#FastQC (Andrews, 2010)  and Trimmomatic (Bolger et al., 2014) 
bsub206 python ../../../scripts/sra.to.fastqc.py -i ./ -o ./ -l ./layout.txt

#for Cvi
java -jar Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ERR754055_1.fastq.gz ERR754055_2.fastq.gz ERR754055_paired_1.fastq.gz ERR754055.unpaired1.fq.gz ERR754055_paired_2.fastq.gz ERR754055_unpaired2.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36

java/jdk1.8.0_60/bin/java -jar Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 SRR34602Merge.fastq.gz SRR34602Merge_trim.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36 &


################################################################### workflow for protein-coding gene annotation ########################################################################
## step 1: run pipeline for protein-coding gene annotation
  bsuball1 python ../scripts/evm.pasa.integrate.pipeline.py -f ./annotation.config
  Read the file "/AMPRILdenovo/annotation/Cvi/run.log" to check the details of each step
	#1) protein sequence alignment using exonerate, protein data are
	exonerate --ryo 'AveragePercentIdentity: %pi\n' --showvulgar no --showalignment no --showquerygff no --model protein2genome --showtargetgff yes  --percent 70 --minintron 10 --maxintron 60000 /AMPRILdenovo/annotation/Cvi/protein/exonerate/prot.chunk.1.fa /AMPRILdenovo/annotation/Cvi/reference/chr.all.v1.0.fasta > /AMPRILdenovo/annotation/Cvi/protein/exonerate/ex.chunk.1.out
	#2) RNAseq reads mapping using hisat2
	/AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/run.ERR754075.sh
	hisat2 --dta --summary-file /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.summary -p 10 -x /AMPRILdenovo/annotation/Cvi/reference/chr.all.v1.0.fasta -1 /AMPRILdenovo/annotation/data/RNAseq/Cvi/ERR754075_paired_1.fastq.gz -2 /AMPRILdenovo/annotation/data/RNAseq/Cvi/ERR754075_paired_2.fastq.gz -S /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.sam
	samtools view -bS /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.sam -o /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.bam
	samtools sort -@ 10  -o /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.srt.bam /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.bam
	stringtie /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.srt.bam -p 10 -o /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.gtf
	echo /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.gtf >> /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/mergelist.txt
	rm /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.sam /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/ERR754075.bam
	
	#3) ab initio prediction using AUGUSTUS, GlimmerHMM, SNAP
	augustus --UTR=on --print_utr=on --exonnames=on --codingseq=on --genemodel=complete --alternatives-from-evidence=true --gff3=on --species=arabidopsis --extrinsicCfgFile=/AMPRILdenovo/annotation/Cvi/abinitio/extrinsic.cfg --outfile=/AMPRILdenovo/annotation/Cvi/abinitio/augustus/augustus.no-hint.chunk.1.gff /AMPRILdenovo/annotation/Cvi/reference/split/genome.chunk.1.fa
	glimmerhmm /AMPRILdenovo/annotation/Cvi/reference/split/genome.chunk.1.fa -d /home/jiao/softwares/genepred/GlimmerHMM/trained_dir/arabidopsis -o /AMPRILdenovo/annotation/Cvi/abinitio/glimmer/glimmer.chunk.1.gff -g -f
	snap /bin/SNAP/snap-2013-11-29/HMM/A.thaliana.hmm /AMPRILdenovo/annotation/Cvi/reference/split/genome.chunk.1.fa -gff > /AMPRILdenovo/annotation/Cvi/abinitio/snap/SNAP.chunk.1.gff
	
	#4) merge result of each kind of evidence 
exonerate2hints.pl --in=/AMPRILdenovo/annotation/Cvi/protein/exonerate.out --source=P --out=/AMPRILdenovo/annotation/Cvi/protein/exonerate.hints --minintronlen=10 --maxintronlen=50000
/bin/EVM_r2012-06-25/EvmUtils/misc/exonerate_gff_to_alignment_gff3.pl /AMPRILdenovo/annotation/Cvi/protein/exonerate.out > /AMPRILdenovo/annotation/Cvi/protein/exonerate.4evm.gff
/bin/TransDecoder/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/stringtie.merged.gtf > /AMPRILdenovo/annotation/Cvi/RNAseq/hisat2/rnaseq.4evm.gff
/bin/EVM_r2012-06-25/EvmUtils/misc/augustus_to_GFF3.pl /AMPRILdenovo/annotation/Cvi/abinitio/augustus.merge.out > /AMPRILdenovo/annotation/Cvi/abinitio/augustus.4evm.gff
/bin/EVM_r2012-06-25/EvmUtils/misc/glimmerHMM_to_GFF3.pl /AMPRILdenovo/annotation/Cvi/abinitio/glimmerhmm.merge.out > /AMPRILdenovo/annotation/Cvi/abinitio/glimmerhmm.4evm.gff
/bin/EVM_r2012-06-25/EvmUtils/misc/SNAP_to_GFF3.jiao.pl /AMPRILdenovo/annotation/Cvi/abinitio/snap.merge.out > /AMPRILdenovo/annotation/Cvi/abinitio/SNAP.4evm.gff
/bin/EVM_r2012-06-25/EvmUtils/misc/SNAP_to_GFF3.jiao.pl /AMPRILdenovo/annotation/Cvi/abinitio/snap.merge.out > /AMPRILdenovo/annotation/Cvi/abinitio/SNAP.4evm.gff

	#5) run EVM + PASA 
	/bin/EVM_r2012-06-25/EvmUtils/.././evidence_modeler.pl -G chr.all.v1.0.fasta -g abinitio.4evm.gff -w /annotation/Cvi/EVM_PASA/evm.weights.txt -e rnaseq.4evm.gff  -p exonerate.4evm.gff  --exec_dir /AMPRILdenovo/annotation/
Cvi/EVM_PASA/split/chr1/chr1_1-100000 > /AMPRILdenovo/annotation/Cvi/EVM_PASA/split/chr1/chr1_1-100000/evm.out


for k in {gene,CDS,prot};do /bin/EVM_r2012-06-25/EvmUtils/gff3_file_to_proteins.pl evm.all.gff3 ../reference/chr.all.v2.0.fasta $k >evm.annotation.$k.fasta ;done &


#### step 2: repeat annotation
	# /AMPRILdenovo/annotation/Cvi/repeat/

  RepeatMasker -species arabidopsis -gff -dir ./ -pa 20 ../../reference/chr.all.v1.0.fasta
  perl ../../../scripts/repeat.classfied.gff3.pl ./chr.all.v1.0.fasta.out.gff ./chr.all.v1.0.fasta.out ./chr.all.v1.0.fasta.repeats.ann.gff3 repeat.ann.stats &
  egrep -v 'Low|Simple|RNA|other|Satellite' chr.all.v1.0.fasta.repeats.ann.gff3 |cut -f 1,4,5,9 >chr.all.v1.0.TE.bed


  perl ../../../scripts/repeat.classfied.gff3.pl ./chr.all.v2.0.fasta.out.gff ./chr.all.v2.0.fasta.out ./chr.all.v2.0.fasta.repeats.ann.gff3 repeat.ann.stats &
  egrep -v 'Low|Simple|RNA|other|Satellite' chr.all.v2.0.fasta.repeats.ann.gff3 |cut -f 1,4,5,9 >chr.all.v2.0.TE.bed

#### step 3: identify TE-related genes
  nohup perl ../../../scripts/remove.TErelated.genes.pl ../../EVM_PASA/evm.annotation.protein.fasta ../../EVM_PASA/evm.annotation.gene.fasta ../RepeatMasker/chr.all.v1.0.TE.bed ../../EVM_PASA/evm.all.gff3 ./ > remove.TE.log & 
  nohup perl ../../scripts/gene.model.gff.stats.pl ./annotation.genes.gff ../reference/chr.all.v1.0.fasta ./gene.model.stats >np.log&

#### step 4: non-coding gene annotation
  nohup cmscan --cpu 20 --tblout Rfam.scan.out ../../../data/Rfam/Rfam.cm ../chr.all.v1.0.split.fasta > run.log  (wd: /AMPRILdenovo/annotation/Eri/noncoding/InfernalRfam
  nohup perl  ../../../scripts/noncoding.infernal.output.parser.pl ./Rfam.scan.out ./ >nohup.out &
  nohup perl  ../../scripts/get.araport11.lncRNA.seq.pl ../../../tair10/TAIR10_chr_all.new.id.fas ../../../tair10/Araport11/Araport11_GFF3_genes_transposons.201606.gff ../../../tair10/Araport11/ >np.log &

#### step 5: get the initial version
## nomenclature follow the TAIR style: https://www.arabidopsis.org/portals/nomenclature/guidelines.jsp
   nohup python ./scripts/gene.id.update.py -i ./ -v 1.0 >log/gene.id.update.log &


#### step 6: evaluation
  awk '{if ($3!="miRNA_primary_transcript" && $3!="pseudogenic_exon" && $3!="pseudogenic_transcript" && $3!="pseudogenic_tRNA" && $3!="transposon_fragment" && $3!="mRNA" && $3!="protein" && $3 !="CDS" && $3!="exon" && $3!="five_prime_UTR" && $3!="three_prime_UTR" && $3!="lnc_RNA") print}' Araport11_GFF3_genes_transposons.201606.gff |grep -v '^#' |cut -f 1,4,5,9 |grep -vP '\.\d+\;Parent' |grep -v  'ChrC|ChrM' >Araport11_gene.TE.chr1-5.bed &

cd /AMPRILdenovo/annotation/Sha/evaluation
mkdir blastnCol misannotation update 

#1) gene seq. blastn
cd blastnCol
ln -s ../../version/Sha.v1.0.gene.fasta ./gene.fasta
blastn -query ./gene.fasta  -db ../../../../tair10/TAIR10_chr_all.fas -num_threads 20 -evalue 1e-5 -out gene.blastn.Col.out &
#nohup blastn -query ../../repeat/TErelated/evm.annotation.gene.fasta -db ../../../../tair10/TAIR10_chr_all.new.id.fas -num_threads 10 -evalue 1e-5 -out gene.blastn.Col.out &

perl annotation.evaluation.by.geneBlastnCol.pl ./gene.blastn.Col.out
awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' ../blastnCol/gene.blastn.besthit.out |sed 's/chloroplast/ChrC/g' |sed 's/mitochondria/ChrM/g' |awk '{if (!/Chr/) print "Chr"$0;else print}' |sed 's/ChrNone/None/g' >query.gene.blastn.besthit.out

#v2.0
nohup blastn -query ../../../../tair10/Araport11/Araport11.prot.genomic.seq.fasta -db ../../version/chr.all.v2.0.fasta -out arap11.prot-gene.blastn.out -evalue 1e-5 &
nohup perl ../../../scripts/assembly.eval.arap11.single.pl ./arap11.prot-gene.blastn.out ../../../../tair10/Araport11/Araport11.prot.genomic.seq.fasta ./araport11.gene.blastn.assembly.out &

#2) prot seq. blastp and gene family clustering
cd /AMPRILdenovo/genefamily/blastpAraport11/Sha
mkdir prot orthomcl; awk '{if (/>/) print $1;else print}' ../../../annotation/Sha/version/Sha.1.0.protein.fasta  |sed 's/>/>Sha|/g' > prot/Sha.fasta
cd prot/; ln -s ../../Col.fasta ./ ; cd ../ ; cat prot/*.fasta > Col-Sha.fasta ; makeblastdb -in ./Col-Sha.fasta -dbtype prot ;sed 's/Kyo/Sha/g' ../Kyo/orthomcl/orthomcl.config  >orthomcl/orthomcl.config; 
blastp -query Col-Sha.fasta -db Col-Sha.fasta -num_threads 40 -evalue 1e-10 -outfmt 6 -out blastout ; orthomclInstallSchema ./orthomcl/orthomcl.config install.sql.log; grep -P "^[^#]" blastout > blastresult; orthomclBlastParser blastresult prot > orthomcl/similarSequences.txt; perl -p -i -e 's/\t(\w+)(\|.*)orthomcl/\t$1$2$1/' orthomcl/similarSequences.txt; perl -p -i -e 's/0\t0/1\t-181/' orthomcl/similarSequences.txt; cd orthomcl ; orthomclLoadBlast ./orthomcl.config similarSequences.txt ; orthomclPairs ./orthomcl.config orthomcl_pairs.log cleanup=all  ; orthomclPairs ./orthomcl.config orthomcl_pairs.log cleanup=no; orthomclDumpPairsFiles ./orthomcl.config ;  mcl mclInput --abc -I 1.5 -o mclOutput -te 20; orthomclMclToGroups group 1 < mclOutput > groups.txt &


#3) find mis-merging, mis-spliting, missing, mis-annotated genes
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

## input files: 
#	Araport11 and accession protein region bed files and sequences files.
       grep gene  ../../repeat/TErelated/annotation.genes.gff|cut -f 1,4,5,9 |sed 's/TU/model/g' >query.prot.gene.bed
     for k in {C24,Cvi,Eri,Kyo,Ler};do grep gene  $k/repeat/TErelated/annotation.genes.gff|cut -f 1,4,5,9 |sed 's/TU/model/g' >$k/evaluation/misannotation/query.prot.gene.bed;done
#	Blastp result of accession proteins against Araport11 proteins
#	OrthoMCL clustering result between accession and Araport11 proteins
         grep AT  ../../../../genefamily/AssV2tmp/An-1/prot/Results_Aug11/Orthogroups.txt |grep evm >groups.txt
#	Blastn result of Araport11 gene sequences against the accession assembly (Blastn result 1)
         awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' Col.prot.besthit.out |grep -P 'AT\d' >Col.prot.besthit.out2
#	Blastn result of accession gene sequences against Col-0 genome sequences. (Blastn result 2)
           awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' query.prot.besthit.out >query.prot.besthit.out2
## output files:
#		potential.mis-merged.gene.txt (function: findMisMer, based on the result of blastp between Col and the Accession)
#		potential.mis-spliting.gene.txt (function: findMisSplit. based on the result of blastp between Col and the Accession)	
#		potential.query.un-assembled.gene.txt (findMissingMisann. based on the result of Col gene blastn against the accession's assembly)
#		potential.missing.gene.txt (findMissingMisann. based on the result of Col gene blastn against the accession's assembly)
#		potential.mis-exon-intron.gene.txt
#		potential.mis-split.gene.by.blastn.txt
#		potential.mis-merge.gene.by.blastn.txt
#		potential.m-vs-m.toBeChecked.by.blastn.txt
#		futher.check.list
#		potential.mis-annotated.gene.txt
#		Araport11.ungrouped.gene.analysis.stat
#		Araport11.ungrouped.gene.analysis.txt
#		Acc(xxx).ungrouped.gene.analysis.stat
#		Acc(xxx).ungrouped.gene.analysis.txt
#		xxx.genes.to.be.updated.added.txt (based on the potential.xxx.txt from blastn-based analysis)
#		xxx.genes.to.be.updated.added.srt.txt	
#
##some cutoffs
    minId = 80
    maxIdf = 10
    maxSplitCov = 0.8
## to find genes:
#	mis-merged (Blastp result, blastn)
#	mis-split (Blastp result, blastn)
#	wrong exon-intron structure (blastn)
#	false protein-coding genes (not annotated in Araport 11 but actually they were assembled in the Col-0 genome) (Blastn result 2)
#	missing genes (not annotated in accession, but actually they were assembled)  (Blastn result 1)

   python ../../../scripts/annotation.evaluate.find-mis.py -g groups.txt -o ./ -r Col -a Sha -n Araport11.gene.blastn.besthit.out -p blastresult -s Araport11.protein.bed -q query.protein.bed -x Araport11.protein.fasta -y query.protein.fasta -c query.gene.blastn.besthit.out &

   nohup python -u ../../../scripts/annotation.evaluate.find-mis.py \
   -g ./groups.txt \
   -o ./ \
   -n  Col.prot.besthit.out2 \
   -c query.prot.besthit.out2 \
   -p blastp.result \
   -s Col.prot.gene.bed \
   -q query.prot.gene.bed \
   -x Col.prot.fasta \
   -y query.prot.fasta  >nohup.out&

 nohup python -u ../../../scripts/annotation.evaluate.find-mis.py -g ./groups.txt -o ./ -n  Col.prot.besthit.out2 -c query.prot.besthit.out2 -p blastp.result -s Col.prot.gene.bed -q query.prot.gene.bed -x Col.prot.fasta -y query.prot.fasta  >nohup.out&
nohup python -u ../../../scripts/annotation.evaluate.find-mis.py -g ./groups.txt -o ./run2 -n  Col.prot.besthit.out2 -c query.prot.besthit.out2 -p blastp.result -s Col.prot.gene.bed -q query.prot.gene.bed -x Col.prot.fasta -y query.prot.fasta -a Col.gene.LoF.txt -b query.gene.LoF.txt -r ../../RNAseq/hisat2/rnaseq.4evm.gff >np.run2&

##################################	step 7: update    #############################################################
##	1) run scipio
  nohup python ../../../scripts/run.scipio.py -i ../../../data/protein/Araport11-split -o ./  -r ../../reference/chr.all.v1.0.fasta >run.log &

  ln -s ../misannotation/Cvi.genes.to.be.updated.added.srt.txt ./; grep chr Cvi.genes.to.be.updated.added.srt.txt >genes.to.be.updated.txt;  grep -v chr Cvi.genes.to.be.updated.added.srt.txt >>genes.to.be.updated.txt ; ln -s ../../abinitio/abinitio.4evm.gff ; cp ../../../An-1/evaluation/update2/srt.gff.pl ./; perl ./srt.gff.pl ./abinitio.4evm.gff SNAP SNAP.ann.gff; perl srt.gff.pl ./abinitio.4evm.gff GlimmerHMM GlimmerHMM.ann.gff  ; cat ../../augustus/genome.chunk.*.gff |egrep -v '#|intron|transcription|codon' > ../../augustus/augustus.ann.gff; ln -s ../../augustus/augustus.ann.gff; ln -s ../misannotation/Araport11.gene.blastn.besthit.out2 ./ ; ln -s ../../../../tair10/Araport11/Araport11_genes.201606.pep.repr.fasta; ln -s ../misannotation/Araport11.protein.bed ./;  ln -s ../../version/Cvi.protein-coding.genes.v1.0.gff ./ ; cat ../../scipio/run2/splitOut/protein.chunk.*.gff >../../scipio/run2/splitOut/scipio.gff; ln -s  ../../scipio/run2/splitOut/scipio.gff ./ ;  ln -s ../../../../wga/results/Cvi/Cvi.wga.snp.indel.gene.LOF.txt ./wga.snp.indel.gene.LoF.txt

##	2) update
#	workdir:
#	input files:
#		gene.to.be.updated.txt
#		xx.protein-coding.genes.v1.0.gff
#		scipio.gff
#		wga.snp.indel.gene.LoF.txt
#		augustus.ann.gff; SNPA.ann.gff; GlimmerHMM.ann.gff
#		Col gene blastn best hit out
#		Col Araport11 protein sequence and bed files
#		ChrCM.txt (a few of organella contigs were not removed in the assembly process)
#	output files:
#		genes.to.be.updated.txt2
#		updated.gff
#		updated.rmdup.gff
#		updated.highConf.gff
#		updated.highConf.prot.fasta
#	Method:
#		1) check the LoF information resulting from WGA-based SNPs and InDel annotation, and the Col protein sequence alignment result from Scipio, add the update information :
#			ChrCM: ChrC or ChrM genes
#			LowConf: low confident genes
#			unchange: keep the previous annotation, 
#			ChangeSci: annotate based on Scipio result (checkScipio: start codon, stop codon, splice sites, frame-shift, premature stop-codon gain; check AugGenes snapGenes glimGene, check AugGenes2[ab initio], check GeneWise)
#			ChangeAug: annotate based on Augustus ( check AugGenes snapGenes glimGene, check AugGenes2[ab initio], check GeneWise)
#			ChangeSciAug: annotate based on scipio and augustus (checkScipio, check AugGenes snapGenes glimGene, check AugGenes2[ab initio], check GeneWise)
#			not-add: not annotate	
#			
#		2) prepare the other annotation result from Augustus-evidence-based method, SNAP ab initio and GlimmerHMM ab initio result
	 cat ../../abinitio/augustus/augustus.hint.chunk.*.gff |egrep -v '#|intron|transcription|codon' > ../../augustus/augustus.ann.gff
#		3) update 
  nohup python ../../../scripts/update.misann.genes.py -u genes.to.be.updated.txt -g An-1.protein-coding.genes.v1.0.gff -o ./ -s scipio.gff -w wga.snp.indel.gene.LoF.txt -c ChrCM.txt -a augustus.ann.gff -n SNAP.ann.gff -l GlimmerHMM.ann.gff -b ../misann2/Araport11.gene.blastn.besthit.out2  -f ../../reference/chr.all.v1.0.fasta -p  Araport11_genes.201606.pep.repr.fasta -i ./Araport11.protein.bed >update.log2 &
  python ../../scripts/annotation.gene.ID.update.py -i update2/updated.highConf.gff -n ../version/An-1.genes.annotation.v1.0.gff -o ../version -v v1.5 -a An-1 -g ../reference/chr.all.v1.0.fasta &

  nohup python -u ../../../scripts/update.misann.genes.py -u genes.to.be.updated.txt -g annotation.genes.gff -o ./ -s scipio.gff -w wga.snp.indel.gene.LoF.txt -c ChrCM.txt -a augustus.ann.gff -n SNAP.ann.gff -l GlimmerHMM.ann.gff -b ./Col.gene.blastn.besthit.bed  -f ../../reference/chr.all.v2.0.fasta -p Col.prot.fasta -i ./Col.prot.gene.bed >update.log &

nohup python -u ../../../scripts/update.misann.genes.py -u genes.to.be.updated.txt -g annotation.genes.gff -o ./run2 -s scipio.gff -x Col.gene.LoF.txt -y query.gene.LoF.txt -c ChrCM.txt -a augustus.ann.gff -n SNAP.4evm.gff -l glimmerhmm.4evm.gff -b ./Col.gene.blastn.besthit.bed  -f ../../reference/chr.all.v2.0.fasta -p Col.prot.fasta -i ./Col.prot.gene.bed >update.run2.log  &
############################## step 8: re-evaluate based on orthomcl gene familiy analysis ######################################################
  grep -v Low updated.rmdup.gff >updated.highConf.gff
  nohup  EVM_r2012-06-25//EvmUtils/gff3_file_to_proteins.pl updated.rmdup.gff ../../reference/chr.all.v1.0.fasta prot >updated.prot.fasta

  perl temp.pl updated.highConf.prot.fasta >bad.prot.log

  perl EVM_r2012-06-25//EvmUtils/gff3_file_to_proteins.pl updated.highConf.gff ../../reference/chr.all.v1.0.fasta prot > updated.highConf.prot.fasta

  ln -s ../../../../../../genefamily/blastpAraport11/C24/prot/Col.fasta ./
  awk '{if (/>/) print $1;else print}' ../updated.highConf.prot.fasta |sed 's/>/>An1|/g' >prot/An1.fasta
  cat prot/*.fasta >Col.C24.fasta
  makeblastdb -in Col.C24.fasta -dbtype prot &
  nohup blastp -query Col.C24.fasta -db Col.C24.fasta -num_threads 60 -evalue 1e-10 -outfmt 6 -out blastresult &
  blastp -query Col.An1.fasta -db Col.An1.fasta -out orthomcl.blastout -evalue 1e-10 -outfmt 6 -num_threads 40  (wd: /AMPRILdenovo/annotation/An-1/evaluation/update2/orthomcl


###gene model stats
 for k in {An-1,C24,Cvi,Eri};do  python -u ./scripts/annotation.gene.ID.update.py -i  $k/evaluation/update/run2/updated.highConf.gff -n ./$k/version2/tmp/$k.genes.annotation.2.0.tmp.gff -o ./$k/version2 -v v2.0 -a $k -g $k/reference/chr.all.v2.0.fasta > $k/version2/updateID.log & done &





###################################################### centromere, telomere, rDNA cluster  ########################################################### 
#workdir : /AMPRILdenovo/annotation/xxx/repeat/trf

#centromere
  nohup bash ./centromere.find.sh > np.log
	perl ../../../scripts/trf.result.parser.pl ./chr.all.v1.0.fasta.1.1.2.80.5.200.2000.dat 80 5 ./ >np.log
	perl  ../../../scripts/centromere.find.pl ./trf.flt.rm.out2 200 90 0.9 ./ >nohup.out 
	awk '{if ($6>50) print ">"$1"_"$2"_"$3"_"$6"\n"$10}' trf.flt.rm.out2 >trf.flt2.unit50.fasta
	##
	ln -s cluster.1.fasta centromere.unit.candidate.fasta #manually check

	makeblastdb -in trf.flt2.unit50.fasta -dbtype nucl 
	blastn -query ./centromere.unit.candidate.fasta -db ./trf.flt2.unit50.fasta -outfmt 6 -out trf.flt2.unit50.blast.cent.candidate.m6 -num_threads 10

	perl  ../../../scripts/centromere.position.predict.pl ./trf.flt2.unit50.blast.cent.candidate.m6 90  ./trf.flt2.unit50.blast.cent.candidate.position
	awk '{if ($2~/chr/ ) print}' ./trf.flt2.unit50.blast.cent.candidate.position |sort -k2,2 -k3,3n |cut -f 2-4|sort -u >centromere.repeat.bed

#telomere
awk '{if ($4==7) print}' trf.flt.rm.out2 >telomere.out
for k in {C24,Cvi,Eri,Kyo,Ler,Sha};do awk '{if ($4==7) print}' $k/repeat/trf/trf.flt.rm.out2 >$k/repeat/trf/telomere.out ;done

#rDNA cluster
#/AMPRILdenovo/annotation/An-1/noncoding/InfernalRfam
awk '{if ($3=="rRNA") print}' noncoding.gene.gff3 >rRNA.gff3
awk '{if (NR<6) print "chr"$0;else print}' ../../version/chr.all.length.v1.0.txt >ass.leng.txt
perl ../../../scripts/rRNA.cluster.pl ./ass.leng.txt ./rRNA.gff3 ./ 10000 1000 >np.log &

for k in {C24,Cvi,Eri,Kyo,Ler,Sha};do perl ./scripts/rRNA.cluster.pl ../assembly/$k/version/ass.v2.leng.txt ./$k/cytogenetic/rDNA/rRNA.gff3 ./$k/cytogenetic/rDNA/ 10000 1000 ../assembly/$k/version/chr.all.v2.0.gap.bed >$k/cytogenetic/rDNA/np.log ; done









