#!/usr/bin/env python
# encoding: utf-8

'''
Created on Sep 16, 2017
@author: jiao@mpip.mpg.de

## splitting ERROR  (two genes blastp with one Col-0 gene, both has high identity)
## merging ERROR (one gene blastp with two Col-0 gene, both has high identity)
## missing genes (good blastn hit, and protein exonerate hit or RNA-seq support)
## wrong gene structure (if scipo gene model has ortholog??)

'''

import sys
import os
import getopt
import re
from Bio import SeqIO

def main(argv):
    groupFile = ""
    outdir = ""
    ref = ""
    alt = ""
    blastnRes = ""
    blastpRes = ""
    refBed = ""
    altBed = ""
    refProtFas = ""
    altProtFas = ""
    #ara11bed = ""
    blastnResCol = ""
    LofRefFile = ""
    LofAltFile = ""
    rnaGff = ""
    try:
        opts, args = getopt.getopt(argv,"g:o:n:p:s:q:x:y:c:a:b:r:",["group=","outdir=","blastn=","blastp=","refBed=","altBed=","refProt=","altProt=","blastnCol=","LofA=","LoFB=", "rna="]) 
    except getopt.GetoptError:
        print 'annotation.evaluate.find-mis.py -g <group> -o <outdir> -n <blastn> -p <blastp> -s <refBed> -q <altBed> -x <refProt> -y <altProt>  -c <blastnCol> -a <LofA> -b <LofB> -r <rna>'
        sys.exit(2)
    if len(opts) == 0 :
        print 'annotation.evaluate.find-mis.py -g <group> -o <outdir> -n <blastn> -p <blastp>  -s <refBed> -q <altBed> -x <refProt> -y <altProt> -c <blastnCol> -a <LofA> -b <LofB> -r <rna>'
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print 'annotation.evaluate.find-mis.py -g <group> -o <outdir>  -n <blastn> -p <blastp> -s <refBed> -q <altBed> -x <refProt> -y <altProt>  -c <blastnCol> -a <LofA> -b <LofB> -r <rna>'
            sys.exit()
        elif opt in ("-g", "--group"):
            groupFile = arg           
        elif opt in ("-o", "--out"):
            outdir = arg
        elif opt in ("-n", "--blastn"):
            blastnRes = arg
        elif opt in ("-p", "--blastp"):
            blastpRes = arg
        elif opt in ("-s", "--refBed"):
            refBed = arg
        elif opt in ("-q", "--altBed"):
            altBed = arg
        elif opt in ("-x", "--refProt"):
            refProtFas = arg
        elif opt in ("-y", "--altProt"):
            altProtFas = arg
        elif opt in ("-c"," --blastnCol") :
            blastnResCol = arg
        elif opt in ("-a"," --LofA") :
            LofRefFile = arg     
        elif opt in ("-b"," --LofB") :
            LofAltFile = arg
        elif opt in ("-r","--rna") :
            rnaGff = arg
    ##some cutoffs
    minId = 80
    maxIdf = 10
    maxSplitCov = 0.8
    
    if not os.path.isdir(outdir) :
        os.system("mkdir -p " + outdir)
    ##
    print refProtFas
    refLen = getProtLeng(refProtFas, 1)
    altLen = getProtLeng(altProtFas, 2)
    
    [refIdx, refPos]= getGeneIndex(refBed,1)
    [altIdx, altPos]= getGeneIndex(altBed,0) 
    
    refGeneBed = getBed(refBed)
    altGeneBed = getBed(altBed)
    
    refBest = getBesthit(blastnRes,1)
    altBest = getBesthit(blastnResCol,0)
    
    ## get specific/un-grouped genes 
    [ungrRef,ungrAlt,ort ] = getSpecGene(groupFile,refLen,altLen)
            
    '''
      Method 1: evaluation based on blastp [Col-Acc] result and Col-gene-acc_assembly blastn result      
    '''
    #check mis-merged genes
    blastpFlt1 = outdir + "/blastp.flt.out1"
    fltBlastp(blastpRes, blastpFlt1, refIdx, altIdx, altLen, refLen)
    [misMergeGene,misMerGeneRef,misMerGeneAlt]= findMergeBlastp(blastpFlt1,outdir,minId,maxSplitCov,maxIdf,ungrRef,ungrAlt,refIdx,altIdx,altLen,altPos, refLen)
    
    #check mis-splitting genes
    [misSplitGene,misSplitGeneRef,misSplitGeneAlt]= findSplitBlastp(blastpFlt1,outdir,minId,maxSplitCov,maxIdf,ungrRef,ungrAlt,refIdx,altIdx,refLen,altPos)
    
    #check unassembled genes ; un-annotated genes; mis-annotated genes (wrong exon-intron structure)   
    print LofRefFile
    colLoF = getColLof(LofRefFile)
    accLoF = getAccLof(LofAltFile)
    
    
    '''
      Method 2: evaluation based on  Col-gene-acc_assembly blastn result          
    1. genes: mis-merged 
    2. genes: mis-split 
    3. genes: mis-exon-intron stucture
    4. genes: mis-annotated (not annotated in Col) 
    5. genes: not annotated in accession assembly, to be added     
    '''    
    
    blastpFlt2 = outdir + "/blastp.flt.out2" 
    [blastp, bestBlastp] = getBlastp(blastpRes, blastpFlt2, refLen, altLen)
    
    [missingRef,misAnnRef,misAnnAlt,missAssRef] = findMM(blastnRes,altBed,outdir,ort,ungrRef,ungrAlt,misMergeGene,misSplitGene,refIdx,altIdx, blastp)
    print "\nPotential missing genes in accession but annotated in Col, resulted ungrouped Col ", len(missingRef)
    k0 = 0
    for g in missingRef:
        if colLoF.has_key(g) :
            k0 += 1
    print "\nPotential missing genes in accession but annotated in Col, resulted ungrouped Col [LoF in other acc]", k0
    
    outfile = outdir + "/potential.Qry.mis-annotated.gene.txt"
    getAltSpeGene(ungrAlt,altPos,outfile)  ## potential specific gene or mis-annotated prot-gene           
    
    [missAssAlt,missingAlt,gtype, rnaSup] = findColMissGene(outdir,blastnResCol, refBed, altBed, ungrAlt, ungrRef, accLoF, rnaGff, blastp)
    
    
    inFile = outdir + "/potential.Qry.mis-merged.gene.by.blastn.txt"
    [misMer, misMerRef] = findMerBlastn(inFile, ungrRef, ungrAlt)
    
    inFile = outdir + "/potential.Qry.mis-split.gene.by.blastn.txt"
    [misSplit, misSplitRef] = findSplitBlastn(inFile, ungrRef, ungrAlt)
    
    inFile = outdir + "/potential.Qry.mis-exon-intron.gene.txt"
    [misExon, misExRef,misExAlt, missingGrAlt]  = findExBlastn(inFile, ungrRef, ungrAlt)
    
    accBlastnCol = getAccBlastnCol(blastnResCol)
    grNum = getGroupNum(groupFile)
    



    ''' step 2a: ungrouped Col gene analysis '''

    fi = open(refBed,"r")
    ara11 = {}
    while True :
        line = fi.readline()
        if not line : break
        if not re.search("protein", line) : continue
        t = line.strip().split("\t")
        id = t[3].split(";")[0]
        id = id.split("=")[1]
        ara11[id] = 1
    fi.close()
    
    fo = open(outdir + "/Col.ungrouped.gene.analysis.txt","w")
    fos = open(outdir + "/Col.ungrouped.gene.analysis.stats","w")
    anaTy1 = ["mis-merge-blastp","mis-split-blastp", "mis-merge-blastn", "mis-split-blastn","missing-ann","mis-exon-intron","not-assembled","obsolete","unknown","unknown-Lof","LoF","mis-ungr"]  
    ana1 = [0,0,0,0,0,0,0,0,0,0,0,0]    
    colUnknown = {}
    for gene in sorted(ungrRef.keys()) :
        flag = "xxx"
        lof = "No-LoF"
        if colLoF.has_key(gene) :
            ana1[10] += 1
            lof = colLoF[gene]
        if not ara11.has_key(gene) :
            flag = "obsolete"
            ana1[7] += 1
        elif misMerGeneRef.has_key(gene) :
            if misMerRef.has_key(gene) :
                flag = "mis-merge-blastp/blastn"
                ana1[2] += 1
            else :
                flag = "mis-merge-blastp"
                #ana1[2] += 1
            ana1[0] += 1
        elif misMerRef.has_key(gene) :
            flag = "mis-merge-blastn"
            ana1[2] += 1            
        elif misSplitGeneRef.has_key(gene) :
            if misSplitRef.has_key(gene) :
                flag = "mis-split-blastp/blastn"
                ana1[3] += 1
            else :
                flag = "mis-split-blastp"
            ana1[1] += 1
        elif misSplitRef.has_key(gene) :
            flag = "mis-split-blastn"
            ana1[3] += 1
                        
        elif missingRef.has_key(gene) :
            flag = "missing-ann"
            ana1[4] += 1
            
        elif misAnnRef.has_key(gene) :
            flag = "mis-exon-intron"
            ana1[5] += 1
            
        elif missAssRef.has_key(gene) :
            if missAssRef[gene] == "ungrouped-unass" :                
                flag = "not-assembled"
            else :
                flag = "partial-assembled"
            ana1[6] += 1
        elif bestBlastp.has_key(gene) and bestBlastp[gene][0] > 80 and bestBlastp[gene][1] > 0.80  and bestBlastp[gene][2] > 0.80 :
            flag = "mis-ungr"
            ana1[11] += 1        
        else :
            flag = "unknown"
            colUnknown[gene]=1
            ana1[8] += 1
            if colLoF.has_key(gene) :
                ana1[9] += 1
        bestp = "xx\txx\txx\txx"
        if bestBlastp.has_key(gene) :
            bestp = "\t".join(str(x) for x in bestBlastp[gene])
        fo.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(refGeneBed[gene], gene, flag, refBest[gene], bestp, lof))
    fo.close()
    
    for k in range(len(ana1)) :
        fos.write("{}\t{}\n".format(anaTy1[k],ana1[k]))
    fos.close()
    
    
    
    
    '''step 2b: analysis of ungrouped genes of AMPRIL other accession genome '''     
    anaTy2 = ["mis-merge-blastp","mis-split-blastp","mis-merge-blastn","mis-split-blastn","missing-ann","mis-exon-intron","not-assembled","xxx","unknown","unknown-LoF","LoF","mis-ungr"]
    ana2 =  [0,0,0,0,0,0,0,0,0,0,0,0]   
        
    fo = open(outdir + "/query.ungrouped.gene.analysis.txt","w")            
    fos = open(outdir + "/query.ungrouped.gene.analysis.stats","w")
    
    for gene in sorted(ungrAlt.keys()) :
        flag = "unknown"
        if misMerGeneAlt.has_key(gene) :
            if misMer.has_key(gene) :                
                flag = "mis-merge-blastp/blastn"
                ana2[2] += 1
            else :
                flag = "mis-merge-blastp"
            ana2[0] += 1
        elif misMer.has_key(gene) :
            flag = "mis-merge-blastn"
            ana2[2] += 1
                        
        elif misSplitGeneAlt.has_key(gene) :
            if misSplit.has_key(gene) :
                flag = "mis-split-blastp/blastn"
                ana2[3] += 1
            else :
                flag = "mis-split-blastp"                
            ana2[1] += 1
        elif misSplit.has_key(gene) :
            flag = "mis-split-blastn"
            ana2[3] += 1
            
        elif missingAlt.has_key(gene) :
            flag = "missing-ann"
            ana2[4] += 1
            
        elif misAnnAlt.has_key(gene) :
            flag = "mis-exon-intron"
            ana2[5] += 1
        elif missAssAlt.has_key(gene) :
            #flag = "not-assembled"
            #ana2[6] += 1       
            if missAssAlt[gene] == "ungrouped-unass" :                
                flag = "not-assembled"
            else :
                flag = "partial-assembled"
            ana2[6] += 1 
        elif bestBlastp.has_key(gene) and bestBlastp[gene][0] > 80 and bestBlastp[gene][1] > 0.80  and bestBlastp[gene][2] > 0.80 :
            flag = "mis-ungr"
            ana2[11] += 1
            
        elif gtype.has_key(gene) :
            geneTy = gtype[gene]
            flag = geneTy
            ana2[7] += 1
            if accLoF.has_key(gene) :
                ana2[9]+=1
        else :
            flag = "unknown"
            ana2[8] += 1
            if accLoF.has_key(gene) :
                ana2[9]+=1
        lof = "No-LoF"
        if accLoF.has_key(gene) :
            ana2[10]+=1
            lof = accLoF[gene]    
        bestp = "xx\txx\txx\txx"        
        if bestBlastp.has_key(gene) :
            bestp = "\t".join(str(x) for x in bestBlastp[gene])                
        fo.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(altGeneBed[gene], gene, flag, altBest[gene], bestp, lof))
    fo.close()
    
    for k in range(len(ana2)) :
        fos.write("{}\t{}\n".format(anaTy2[k],ana2[k]))
    fos.close()




    ''' step 3: final stats: '''
    #acc gene Col-sp Acc-sp
    st1 = [0,0,0] 
    st1[0] = len(altLen.keys())
    st1[1] = len(ungrRef.keys())
    st1[2] = len(ungrAlt.keys())
    #s-blastp s-blastn split m-blastp m-blastn merge ex-in acc.unann acc.unass acc.pa.ass col.unass col.par.ass col.unann Col-LoF Acc-LoF
    st2 = [0,0,0, 0,0,0, 0, 0,0,0, 0,0,0, 0,0]
    st2[0] = len(misSplitGene.keys())
    st2[1] = len(misSplit.keys())
    tmp = misSplit
    tmp.update(misSplitGene)
    st2[2] = len(tmp.keys())
    
    st2[3] = len(misMergeGene.keys())
    st2[4] = len(misMer.keys())
    tmp = misMer
    tmp.update(misMergeGene)
    st2[5] = len(tmp.keys())
    
    st2[6] = len(misExon.keys())
    st2[7] = len(missingRef.keys())
    for x in missAssRef.keys() :
        if re.search(r"unass",missAssRef[x]) :
            st2[8] +=1
        else :
            st2[9] +=1            
    for x in missAssAlt.keys() :
        if re.search(r"unass", missAssAlt[x]) :
            st2[10] +=1
        else :
            st2[11] +=1
    st2[12] = len(missingAlt.keys())
    st2[13] = ana1[10]
    st2[14] = ana2[10]
    
    #Col-sp split merge ex-in acc.unann acc.unass acc.pa.ass unk-LoF
    st3 = [0, 0, 0, 0, 0, 0, 0,0 ]
    st3[0] = st1[1]
    for x in misSplitGeneRef.keys() :
        if misSplitGeneRef[x] == "ungrouped" :
            st3[1] += 1
    for x in misSplitRef.keys() :
        if not misSplitGeneRef.has_key(x) and misSplitRef[x] == "ungrouped":
            st3[1] += 1
            
    
    for x in misMerGeneRef.keys() :
        if misMerGeneRef[x] == "ungrouped" :
            st3[2] += 1
    for x in misMerRef.keys() :
        if not misMerGeneRef.has_key(x) and misMerRef[x] == "ungrouped":
            st3[2] += 1        
    st3[3] = misExRef
    for x  in missingRef.keys() :
        if ungrRef.has_key(x) :
            st3[4] += 1
    for x in missAssRef.keys() :
        if missAssRef[x] == "ungrouped-unass" :
            st3[5] += 1
        elif missAssRef[x] == "ungrouped-parass" :
            st3[6] += 1
    st3[7] = ana1[9]
    
    #Acc-sp split merge ex-in col.unann col.unass col.pa.ass unk-LoF ungr-rnaSup
    st4 = [0, 0, 0, 0, 0, 0, 0,0 ,0 ]
    st4[0] = st1[2]  
    for x in misSplitGeneAlt.keys() :
        if misSplitGeneAlt[x] == "ungrouped" :
            st4[1] += 1
    for x in misMerGeneAlt.keys() :
        if misMerGeneAlt[x] == "ungrouped" :
            st4[2] += 1
    st4[3] = misExAlt
    for x in missingAlt.keys() :
        if ungrAlt.has_key(x) :
            st4[4] += 1 
    for x in missAssAlt.keys() :
        if missAssAlt[x] == "ungrouped-unass" :
            st4[5] += 1
        elif missAssAlt[x] == "ungrouped-parass" :
            st4[6] += 1
    st4[7] = ana2[9]
    for x in rnaSup :
        if ungrAlt.has_key(x) :
            st4[8]+=1
    fo = open(outdir + "/annotation.evaluation.stats","w")
    fo.write("all\t{}\t{}\n".format( "\t".join(str(x) for x in(st1)), "\t".join(str(x) for x in(st2)) ))
    fo.write("\n")
    fo.write("Col-sp\t{}\tAcc-sp\t{}\n".format( "\t".join(str(x) for x in(st3)), "\t".join(str(x) for x in(st4)) ))
    
    




    '''step 4: combine the above information to get query gene to be updated or added '''
    fi = open(altBed,"r")
    outfile1 = outdir + "/query.genes.to.be.updated.added.txt" 
    fo = open(outfile1,"w")
    k0 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[3].split(";")[0]
        id = id.split("=")[1]
        blast = ";".join(accBlastnCol[id])
        
        flag = ""
        gr = "grouped"
        if ungrAlt.has_key(id) : 
            gr = "ungrouped"
        else :
            gr = grNum[id]
        if misMer.has_key(id) :
            refGene = misMer[id][0]
            grRef = misMer[id][1]
            for k in range(len(misMer[id])) :
                if k == 0 : continue
                if k % 2 == 0 :
                    refGene = refGene + "," + misMer[id][k]
                    grRef = grRef + "," + misMer[id][k+1]
            flag = "mis-merge" + "\t" + refGene + "\t" + grRef
        elif misMergeGene.has_key(id) :
            refGene = misMergeGene[id][0]
            grRef = misMergeGene[id][1]
            for k in range(len(misMergeGene[id])) :
                if k == 0 : continue
                if k % 2 == 0 :
                    refGene = refGene + "," + misMergeGene[id][k]
                    grRef = grRef + "," + misMergeGene[id][k+1]
            flag = "mis-merge" + "\t" + refGene + "\t" + grRef    
             
        elif misSplit.has_key(id) :
            refGene = misSplit[id][0]
            flag = "mis-split" + "\t" + refGene + "\t" + misSplit[id][1]
            
        elif misSplitGene.has_key(id):
            refGene = misSplitGene[id][0]
            flag = "mis-split" + "\t" + refGene + "\t" + misSplitGene[id][1]
        
        elif missingGrAlt.has_key(id) :
            refId = missingGrAlt[id]
            flag = "mis-ungr\t" + refId + "\t-"
                
        elif misExon.has_key(id) :
            #flag = "mis-ex-int"
            refGene = misExon[id][0]
            flag = "mis-exon" + "\t" + refGene + "\t" + misExon[id][2]
            
        elif missAssAlt.has_key(id) :
            if missAssAlt[id] == "ungrouped-unass" :
                flag = "Col-Not-ass" + "\t-\t-"
            else :
                flag = "Col-Par-ass" + "\t-\t-"
                        
        elif missingAlt.has_key(id) :
            if rnaSup.has_key(id) :
                flag = "ColNo-RNAsup" + "\t" + ",".join(missingAlt[id]) + "\t-"
            else :
                flag = "Col-No-ann" + "\t" + ",".join(missingAlt[id]) + "\t-"
         
        elif gtype.has_key(id) :
            flag = "Col-diff-ann" + "\t" + gtype[id] + "\t-"
                
        elif ungrAlt.has_key(id) :
            if accLoF.has_key(id) :                
                flag = "Col-LoF" + "\t"+ accLoF[id] + "\t-"
            else :
                flag = "Col-No-ann" + "\t-\t-"   
                k0 += 1                 
        else :
            flag = "unchange\t-\t-"        
        
        fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(t[0],t[1],t[2],id,gr,blast,flag))     
    fi.close()
    print "ungrouped acc , col no LoF and no ann, ", k0
    #not annotated in acc but in Col
    missRefUngr = 0
    for k in missingRef :
        [chrom,start,end,prot,cov,iden] = missingRef[k]
        flag = "add"
        if ungrRef.has_key(k) :
            flag = flag + "\t" + k + "\tungrouped"
            missRefUngr += 1
        else  :
            flag = flag + "\t" + k + "\tgrouped"
            
        fo.write("{}\t{}\t{}\t{}\t{}\t{};{}\t{}\n".format(chrom,start,end,k,"**",cov,iden,flag))
    
    
    fi = open(outdir + "/Col.gene.blastn.intersect.Qry.gene.txt","r")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        refId = t[3].split(".")[0]        
        if colUnknown.has_key(refId) :
            flag = "add"
            if ungrRef.has_key(refId) :
                flag = flag + "\t" + refId + "\tungrouped"
            else  :
                flag = flag + "\t" + refId + "\tgrouped"
            fo.write("{}\t{}\t{}\t{}\t{}\t{};{}\t{}\n".format(t[0],t[1],t[2],refId,"**",t[4],t[5],flag))
    fi.close()
    fo.close()
    outfile2 = outdir + "/query.genes.to.be.updated.added.srt.txt"
    os.system("sort -k1,1 -k2,2n " + outfile1  + " > " + outfile2)
            
def getBesthit(inFile, isCol):
    besthit = {}
    fi = open(inFile,"r")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if isCol == 1 :
            id = t[3].split(".")[0]
            besthit[id] = line.strip()
        else :
            besthit[t[3]] = line.strip()
    fi.close()
    return besthit


def getColLof(inFile):
    colLof = {}
    fi = open(inFile,"r")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[5] != "No" :
            colLof[t[0]] = t[5]
    fi.close()
    return colLof

def getAccLof(inFile):
    accLof = {}
    fi = open(inFile,"r")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[5] != "No" :
            accLof[t[1]] = t[5]
    fi.close()
    return accLof

def getBlastp(inFile,outFile, refLen, altLen):
    fi = open(inFile,"r")
    fo = open(outFile,"w")
    blastp = {}
    bestBlastp = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if t[2] < 60 : continue
        t[1]=t[1].split(".")[0]
        cov1 = abs(float(t[7])-float(t[6]))/altLen[t[0]]
        cov2 = abs(float(t[9])-float(t[8]))/refLen[t[1]]
        if cov1 < 0.60 or cov2 < 0.60 : continue
        fo.write("{}\t{}\t{}\n".format(line.strip(), cov1, cov2))
        refId = t[1].split(".")[0]
        blastp[t[0] + refId] = [t[2], cov1, cov2]
        if not bestBlastp.has_key(t[0]) :
            bestBlastp[t[0]] = [refId, float(t[2]), cov1, cov2]
        else :
            if bestBlastp[t[0]][1]*bestBlastp[t[0]][2]*bestBlastp[t[0]][3] < float(t[2])*cov1*cov2 :
                bestBlastp[t[0]] = [refId, float(t[2]), cov1, cov2]
        
        if not bestBlastp.has_key(refId) :
            bestBlastp[refId] = [t[0], float(t[2]), cov1, cov2]
        else :
            if bestBlastp[refId][1]*bestBlastp[refId][2]*bestBlastp[refId][3] < float(t[2])*cov1*cov2 :
                bestBlastp[refId] = [t[0], float(t[2]), cov1, cov2]        
    fi.close()
    return [blastp, bestBlastp]


def getAccBlastnCol (blastnResCol):
    fi = open(blastnResCol,"r")
    accBlastnCol = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        #id = t[3].split(".")[0]
        id = t[3]
        #if re.search("evm",t[3]) :
        #    id = re.sub("model","TU",t[3]) 
        accBlastnCol[id] = [t[0],t[1],t[2],t[4],t[5]]
    fi.close()
    return accBlastnCol

def getGroupNum (groupFile):    
    fi = open(groupFile,"r")
    grNum = {}
    ref = "AT"
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split()
        refIds = []
        altIds = []        
        for k in range(len(t)) :
            if k == 0 : continue
            if re.search(ref,t[k]) :
                id = t[k].split(".")[0]
                refIds.append(id)
            else :            
                id = t[k]
                altIds.append(id)
        for k in altIds :
            grNum[k] = str(len(refIds)) + "-" + str(len(altIds))      
    fi.close()  
    return grNum

    
def findMerBlastn(inFile, ungrRef, ungrAlt): 
    fi = open(inFile,"r")
    misMer={}
    misMerRef = {}
    k1 = 0
    k2 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        #misMer[t[0]]=[]
        altId = t[9].split(";")[0]
        altId = altId.split("=")[1]
        refId = t[3].split(".")[0]
        gr = "grouped"        
        if ungrRef.has_key(refId) : 
            gr = "ungrouped"
            k1 += 1
        misMerRef[refId] = gr
         
        if ungrAlt.has_key(altId) :
            k2 += 1   
                            
        if not misMer.has_key(altId) :
            misMer[altId]=[refId,gr]
        else :
            misMer[altId].append(refId)
            misMer[altId].append(gr)
    fi.close()
    print "\n find mis-merged gene by blastn"
    print "mis-merged genes acc number: ", len(misMer)
    print "mis-merged genes col number: ", len(misMerRef)
    print "mis-merged genes resulted in ungroupd Col ", k1
    print "mis-merged genes resulted in ungroupd Acc ", k2
    return misMer, misMerRef

def findSplitBlastn(inFile, ungrRef, ungrAlt):    
    fi = open(inFile,"r")
    misSplit = {}
    misSplitRef = {}
    k1 = 0
    k2 = 0 
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        altId = t[9].split(";")[0]
        altId = altId.split("=")[1]
        refId = t[3].split(".")[0]
        misSplitRef[refId] = 1
        gr = "grouped"
        if ungrRef.has_key(refId) : 
            gr = "ungrouped"
            k1 += 1
        if ungrAlt.has_key(altId) :
            k2 += 1
        misSplit[altId]=[refId,gr]
        #misSplit[x]=[t[4],xx2[1],t[6]]
    fi.close()
    print "\n # find mis-split gene by blastn "
    print "mis-split gene acc number :" , len(misSplit)
    print "mis-split gene result in ungrouped Col gene :" , k1
    print "mis-split gene result in ungrouped Acc gene :" , k2
    return [misSplit,misSplitRef]
    

def findExBlastn(inFile, ungrRef, ungrAlt):
    fi = open(inFile,"r")
    misExon={}
    missingGrAlt = {}
    k = 0
    k1 = 0
    k2 = 0
    k31 = 0
    k32 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        #misMer[t[0]]=[]
        id = t[9].split(";")[0]        
        id = id.split("=")[1]
        t[11] = re.sub("ref\-|alt\-","",t[11])
        t[12] = re.sub("ref\-|alt\-","",t[12])        
        refId = t[3].split(".")[0]
        if t[-1] == "not-ortholog" :
            misExon[id] = [refId,t[12],t[11]]        
        if t[-1] == "missing-ortholog" :
            missingGrAlt[id] = refId
            if ungrRef.has_key(refId) :
                k31 += 1
            if ungrAlt.has_key(id) :
                k32 += 1
        else :
            if ungrAlt.has_key(id) :
                k2 += 1
            if ungrRef.has_key(refId) :
                k1 += 1
        #print id,t[12],t[11]
        #sys.exit()
    fi.close()
    print "\n # find gene with wrong exon-intron structure based on blastn"
    print "mis-exon gene acc number", len(misExon)
    print "mis-exon genes resulted in ungrouped Col gene", k1
    print "mis-exon genes resulted in ungrouped Acc gene", k2
    
    print "mis-ungrouped Col gene", k31
    print "mis-ungrouped Acc gene", k32
    print "len(missingGrAlt):", len(missingGrAlt)
    misExRef = k1
    misExAlt = k2
    return [misExon, misExRef,misExAlt, missingGrAlt]      
    
def getBed(bedFile):
    fi = open(bedFile,"r")
    geneBed = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[3].split(";")[0]
        id = id.split("=")[1]
        geneBed[id]="\t".join(t[0:3])
    fi.close()
    return geneBed
            


def findColMissGene(outdir,blastnRes,ara11bed,accBed,ungrAlt, ungrRef, accLoF, rnaGff, blastp):
   # missingRef = {}
    #misAnnRef = {}
    #misAnnAlt = {}
    missAssAlt = {}
    missingAlt = {}
    #missingAlt2 = {}
    geneBed = getBed(accBed)
    '''
      awk '{ if (/None|chl|mito/) print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10 ; else print "Chr"$2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}'
      ../blastnCol/gene.blastn.besthit.out >query.prot.besthit.out2
    '''        
    fi = open(blastnRes,"r") 
    fo = open(outdir + "/potential.Col.un-assembled.gene.txt", "w")
    k1 = 0
    k2 = 0
    while True :
        line = fi.readline()
        if not line :break
        t = line.strip().split("\t")
        id = t[3]
        if t[0] == "None" :
            #id = t[3].split(".")[0]
            if ungrAlt.has_key(id) :                         
                missAssAlt[id] = "ungrouped-unass"
                fo.write("{}\t{}\t{}\n".format(geneBed[id],id,"ungrouped-unass"))
                k1 += 1
            else :
                missAssAlt[id] = "grouped-unass"
                fo.write("{}\t{}\t{}\n".format(geneBed[id],id,"grouped-unass"))
            
        elif float(t[5]) > 90 and float(t[4]) < 90 :
            missAssAlt[id] = 2
            if ungrAlt.has_key(id) :                         
                missAssAlt[id] = "ungrouped-parass"
                fo.write("{}\t{}\t{}\n".format(geneBed[id],id,"ungrouped-parass"))
                k2 += 1
            else :
                missAssAlt[id] = "grouped-parass"
                fo.write("{}\t{}\t{}\n".format(geneBed[id],id,"grouped-parass"))
    fi.close()
    fo.close() 
    print "\n###"
    print "Gene in other accession but not assembled or partially assembled in Col genome (blastn): " , len(missAssAlt)
    print "Gene in other accession but not assembled in Col genome, resulted in ungrouped Acc: " , k1
    print "Gene in other accession but partially assembled in Col genome, resulted in ungrouped Acc: " , k2
    
    blastnRes2 = outdir + "/query.prot-gene.blastn.Col.besthit.out2"
    os.system("grep -v None " + blastnRes + " |sort -k1,1 -k2,2n -k3,3n > " + blastnRes2)
    os.system("intersectBed -a " + blastnRes2 + " -b " + ara11bed + " -v |sort -k1,1 -k2,2n -k3,3n > " + outdir + "/potential.Col.missing.gene.txt")
    fi = open(outdir + "/potential.Col.missing.gene.txt", "r" )
    fo = open(outdir + "/potential.Col.missing.gene.txt2", "w" )
    
    k = 0
    k0 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")        
        #id = t[3].split(".")[0]
        id = t[3]
        #if re.search("evm",t[3]) :
        #    id = re.sub("model","TU",t[3]) 
        fo.write("{}\t{}".format(geneBed[id],line.strip()))
        if re.search("ChrC|ChrM",t[0]) :
            fo.write("\t{}\n".format(t[0]))
        elif ungrAlt.has_key(id) :  
            # missingRef[id] = 1          
            fo.write("\tungrouped\n")
            k += 1
            missingAlt[id] = t
            if accLoF.has_key(id) :
                k0 +=1
        else :
            fo.write("\tgrouped\n")
            missingAlt[id] = t
    fi.close()
    fo.close()
    print "\n###"    
    print "Gene annotated in accession and blasted to Col assembly but not annotated in Col :" , len(missingAlt)
    print "Gene annotated in accession and blasted to Col assembly but not annotated in Col, resulted ungrouped alt :" , k
    print "Gene annotated in accession and blasted to Col assembly but not annotated in Col, resulted ungrouped alt (LoF in Col) :" , k0
    os.system("intersectBed -a " + outdir + "/potential.Col.missing.gene.txt2 -b " + rnaGff + " -wo > " + outdir + "/potential.Col.missing.gene.rna.sup.txt")
    fi = open(outdir + "/potential.Col.missing.gene.rna.sup.txt", "r")
    rnaGene = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if not rnaGene.has_key(t[6]) :
            rnaGene[t[6]] = int(t[-1])
        else :
            rnaGene[t[6]] += int(t[-1])
    fi.close()
    
    fi = open(outdir + "/potential.Col.missing.gene.txt2", "r" )
    fo = open(outdir + "/potential.Col.missing.gene.txt3", "w" )
    rnaSup = {}
    k0 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        flag = "xxx"
        if rnaGene.has_key(t[6]) :
            flag = float(rnaGene[t[6]])/(float(t[2]) - float(t[1]))
            if flag > 0.9 :
                rnaSup[t[6]] = 1
                if accLoF.has_key(t[6]) :
                    k0+=1
        if accLoF.has_key(t[6]) :
            flag = str(flag) + "\t" + "Col-LoF"
        else :
            flag = str(flag) + "\t" + "Col-noLoF"
        fo.write("{}\t{}\n".format(line.strip(), flag))
    fi.close()
    fo.close()
    print "potential col missing gene with RNA-sup", len(rnaSup)
    print "potential col missing gene with RNA-sup and Col has LoF", k0
    
    os.system("intersectBed -a " + blastnRes2 + " -b " + ara11bed + " -wo > " + outdir + "/query.gene.blastn.intersect.Col.gene.txt")
    ###if the paired gene is not ortholog, then the gene may have wrong exon-intron structure. Also check the protein exonerate alignment
    fi = open(outdir + "/query.gene.blastn.intersect.Col.gene.txt", "r")
    fo = open(outdir + "/potential.Qry.mis-annotated.gene.txt2", "w")
    prevAltId = ""
    gtype = {}
    k1 = 0
    k2 = 0
    k3 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        #altId = t[3].split(".")[0]
        altId = t[3]
        #if re.search("evm",t[3]) :
        #    altId = re.sub("model","TU",t[3])
        if not re.search("protein", t[9]) : continue
        if not ungrAlt.has_key(altId) : 
            #k1 += 1
            continue
        if  float(t[-1])/(float(t[2]) - float(t[1])) < 0.8 or  float(t[-1])/(float(t[8]) - float(t[7])) < 0.8:
            continue
        
        refId = t[9].split(";")[0]        
        refId = refId.split("=")[1]
        if not ungrRef.has_key(refId) :
            k2 += 1
        if blastp.has_key(t[3] + refId) :            
            k3+=1
            continue
        
        tmp =re.search("locus_type=([\w\_]+)",t[9])        
        ty = tmp.group(0)
        #print ty,altId
        if prevAltId == "" :
            prevAltId = altId
            fo.write(line)
            gtype[altId] = ty
        elif prevAltId != altId :
            prevAltId = altId
            fo.write(line)
            gtype[altId] = ty
        else :
            gtype[altId] = gtype[altId] + ";" + ty                             
    fi.close()
    fo.close()
    print "ungroupd Acc Genes annotated in accession ad blast hit a gene in Col assembly but they are not orthologs: " , len(gtype)
    #print "Genes annotated in accession ad blast hit a gene in Col assembly, they are not orthologs, but the Alt has ortholog: " , k1
    print "ungrouped Acc Genes annotated in accession ad blast hit a gene in Col assembly, they are not orthologs, but the Col has ortholog: " , k2
    print "ungroupded Genes annotated in accession ad blast hit a gene in Col assembly, they are not orthologs, but they have high similarity: " , k3
    return [missAssAlt,missingAlt,gtype, rnaSup]

    
def getGeneIndex (geneBed, ref):

    fi = open(geneBed,"r")
    geneIdx = {}
    n = 0
    genePos = {}
    while True :
        line = fi.readline()
        if not line : break
        #if (not re.search("Note=protein_coding_gene",line) ) and (not re.search("locus_type=protein",line)) : continue
        if ref == 1 and not re.search(r"protein",line) : continue 
        t = line.strip().split("\t")
        #id = t[4].split(";")[0]
        id = t[3].split(";")[0]
        id = id.split("=")[1]
        geneIdx[id]=n
        genePos[id]=[t[0],t[1],t[2]]
        n += 1
    fi.close()
    
    return [geneIdx,genePos]
        
    ### step 1: find mis-merging genes
def fltBlastp(blastpRes,outFile, refIdx,altIdx,altLen, refLen):
    fi = open(blastpRes,"r")    
    fo = open(outFile,"w")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if float(t[2]) < 60 : continue
        aIdx =  altIdx[t[0]]
        rId = t[1].split(".")[0]
        rIdx = refIdx[rId]
        len1 = altLen[t[0]]
        len2 = refLen[rId]
        cov1 = abs((float(t[7]) - float(t[6])))/float(len1)
        cov2 = abs((float(t[9]) - float(t[8])))/float(len2)
        fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line.strip(), aIdx, rIdx, len1, len2, cov1, cov2))
    fi.close()
    fo.close()
    
def findMergeBlastp (blastpFlt,outdir,minId,maxSplitCov,maxIdf,ungrRef,ungrAlt,refIdx,altIdx,altLen,altPos, refLen):
    
    
    
    blastpFltSrt = outdir + "/blastp.flt.srt1.out"
    os.system("sort -k13,13n -k14,14n -k7,7n " + blastpFlt + "  > " + outdir + "/blastp.flt.srt1.out") ## for mis-merging
    
    
    fi = open(blastpFltSrt,"r")
    fo1 = open(outdir + "/potential.Qry.mis-merged.gene.by.blastp.txt","w")
    #fo2 = open(outdir + "/potential.mis-merged.gene.blastp.result.txt","w")
    aIdx = 0 
    rIdx = 0
    [prevId,preCov] = [0,0]
    prevRefGene = ""
    prevAltEnd = 0
    minDist = 10
    misMerGene = {}
    misMerGeneRef = {}
    misMerGeneAlt = {}
    preLine = ""
    k1 =0
    k2 =0
    ref = "AT"
    alt = "evm"
    print "\n### find mis-merging genes in accession based on the blastp results "
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if re.search(alt,t[0]) and re.search(ref,t[1]) :
            if float(t[2]) <= minId : continue
            
            aId = t[0]
            rId = t[1].split(".")[0]            
            
            if aIdx == altIdx[aId] and refIdx[rId] == rIdx + 1 and abs(float(t[2])- prevId) < maxIdf and preCov < maxSplitCov and float(t[3])/altLen[aId] < maxSplitCov and abs(int(t[6]) - prevAltEnd) < minDist:
                [chrom,start,end]=altPos[aId]
                flag1 = "grouped" 
                flag2 = "grouped" 
                misMerGene[aId]= []
                if ungrAlt.has_key(aId) : 
                    flag1 = "ungrouped"
                    k1 += 1
                misMerGeneAlt[aId] = flag1
                if ungrRef.has_key(prevRefGene) : 
                    flag2 = "ungrouped"
                    k2 += 1
                misMerGeneRef[prevRefGene] = flag2
                
                misMerGene[aId].append(prevRefGene)
                misMerGene[aId].append(flag2)
                    
                misMerGene[aId].append(rId)
                if ungrRef.has_key(rId) : 
                    flag2 = flag2 + ",ungrouped"
                    misMerGeneRef[rId] = "ungrouped"                    
                    misMerGene[aId].append("ungrouped")
                    k2 += 1
                else :
                    flag2 = flag2 + ",grouped"         
                    misMerGene[aId].append("grouped")
                    misMerGeneRef[rId] = "grouped"       
                fo1.write("{}\t{}\t{}\t{}\t{},{}\t{}\t{}\n".format(aId,chrom,start,end,prevRefGene,rId,flag1,flag2))
                fo1.write("\t#{}\n\t#{}\n".format(preLine,line.strip()))
                aIdx = altIdx[aId]
                rIdx = refIdx[rId]
                prevId = 0
                preCov =  0
                prevRefGene = "" 
                prevAltEnd = int(t[7])
                preLine = line.strip()
            else :
                aIdx = altIdx[aId]
                rIdx = refIdx[rId]  
                prevId = float(t[2])
                preCov =  float(t[3])/altLen[aId]
                prevRefGene = rId 
                prevAltEnd = int(t[7])     
                preLine = line.strip()
                    
    fi.close()
    fo1.close()
    print "mismerged gene acc : ",len(misMerGene)
    print "mismerged gene col : ",len(misMerGeneRef)
    print "mismerged gene acc number(ungrouped): ", k1
    print "mismerged gene Col number(ungrouped: ", k2
    print "\n"
    return [misMerGene,misMerGeneRef,misMerGeneAlt]
    
    ## step 2: find mis-spliting genes
def findSplitBlastp (blastpFlt,outdir,minId,maxSplitCov,maxIdf,ungrRef,ungrAlt,refIdx,altIdx,refLen,altPos):
    print "### find mis-split genes in accession based on the blastp results "
    os.system("sort -k14,14n -k13,13n -k9,9n " + blastpFlt + "  > " + outdir + "/blastp.flt.srt2.out") ## for mis-spliting
    fi = open(outdir + "/blastp.flt.srt2.out","r")
    fo = open(outdir + "/potential.Qry.mis-split.gene.by.blastp.txt","w")
    aIdx = 0 
    rIdx = 0
    ref = "AT"
    alt = "evm"
    [prevId,preCov] = [0,0]
    preAltGene = ""
    preRefEnd = 0
    misSplitGene = {}
    misSplitGeneRef = {}
    misSplitGeneAlt = {}
    k1 = 0
    k2 = 0
    #k1k2 = 0
    preLine = ""
    preAltCov = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        if re.search(alt,t[0]) and re.search(ref,t[1]) :
            if float(t[2]) < minId : continue
            aId = t[0]
            rId = t[1]
            rId = rId.split(".")[0]
            if aIdx == altIdx[aId] - 1 and refIdx[rId] == rIdx and abs(float(t[2])- prevId) < maxIdf and preCov < maxSplitCov and float(t[3])/refLen[rId] < maxSplitCov and  float(t[-2]) > maxSplitCov and preAltCov > maxSplitCov :
                [chrom,start,end]=[altPos[aId][0],altPos[preAltGene][1],altPos[aId][2]] 
                flag1 = "grouped" 
                
                if ungrAlt.has_key(preAltGene) : 
                    flag1 = "ungrouped"
                    k1 += 1
                misSplitGeneAlt[preAltGene] = flag1
                if ungrAlt.has_key(aId) :
                    k1 += 1 
                    flag1 = flag1 + ",ungrouped"
                    misSplitGeneAlt[aId] = "ungrouped"
                else :
                    flag1 = flag1 + ",grouped"
                    misSplitGeneAlt[aId] = "grouped"
                                        
                flag2 = "grouped"
                if ungrRef.has_key(rId) : 
                    k2 += 1
                    flag2 = "ungrouped"
                misSplitGeneRef[rId] = flag2
                    
                misSplitGene[aId] = [rId,flag2]
                fo.write("{},{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(preAltGene,aId,chrom,start,end,rId,flag1,flag2))
                fo.write("\t#{}\n\t#{}\n".format(preLine,line.strip()))
                #misSplitGene[preAltGene]=1
                #misSplitGene[aId]=1
                aIdx = altIdx[aId]
                rIdx = refIdx[rId]
                prevId = 0
                preCov =  0
                preAltGene = ""
                preAltCov = float(t[-2])
                preLine = line.strip()
            else :
                aIdx = altIdx[aId]
                rIdx = refIdx[rId]      
                prevId = float(t[2])
                preCov =  float(t[3])/refLen[rId] 
                preAltGene = aId
                preAltCov = float(t[-2])
                preLine = line.strip()      
                                
    fi.close()
    fo.close()
    
    print "mis-split gene acc :", len(misSplitGeneAlt)
    print "mis-split gene col :", len(misSplitGeneRef)
    print "mis-split gene acc number(ungrouped): ", k1
    print "mis-split gene Col number(ungrouped: ", k2
    print "\n"
    return [misSplitGene , misSplitGeneRef, misSplitGeneAlt]
    ## step 3 : find missing genes

def findMM (blastnRes,altBed,outdir,ort,ungrRef,ungrAlt,misMer,misSpl,refIdx,altIdx, blastp):
    ##blastnRes includes gene has blastn hit region in accession genome, but no gene was annotated 
    missingRef = {}
    misAnnRef = {}
    misAnnAlt = {}
    missAssRef = {}
    
    print "\n ### find Col genes which are not assembled (could be deleted) in other accession based on Col-gene blastn result "
    
    fi = open(blastnRes,"r")
    fo = open(outdir + "/potential.Qry.un-assembled.gene.txt","w")
    k1 = 0
    k2 = 0
    x = 0
    y = 0 
    while True :
        line = fi.readline()
        if not line :break
        t = line.strip().split("\t")
        if t[0] == "None" :
            #id = t[4].split(".")[0]
            id = t[3].split(".")[0]
            flag = "ungrouped"
            if not ungrRef.has_key(id) : 
                flag = "grouped"
            else :
                k1 += 1
            fo.write("{}\t{}\t{}\t{}\n".format(t[3],"un-assembled",flag, line.strip() ))
            missAssRef[id] = flag + "-unass"
            x += 1
        elif float(t[5]) > 90 and float(t[4]) < 90 :
            y += 1
            id = t[3].split(".")[0]
            flag = "ungrouped"
            if not ungrRef.has_key(id) : 
                flag = "grouped"
            else :
                k2 += 1            
            fo.write("{}\t{}\t{}\t{}\n".format(t[3],"partially-assembled",flag, line.strip() ))
            missAssRef[id] = flag + "-parass"
    fi.close() 
    fo.close()
    print "Unassembled or deleted Col gene in other accession: ", x
    print "Unassembled or deleted Col gene, resulted in Col ungrouped: ", k1
    print "\n"
    print "partially assembled [cov<0.9] or deleted Col gene in other accession: ", y
    print "partially [cov<0.9] unassembled or deleted Col gene, resulted in Col ungrouped: ", k2
    
    blastnRes2 = outdir + "/Col.gene.blastn.besthit.bed"
    os.system("grep -v None " + blastnRes + " > " + blastnRes2)    
    
    os.system("intersectBed -a " + blastnRes2 + " -b " + altBed + " -wao > " + outdir + "/Col.gene.blastn.intersect.Qry.gene.txt")
    os.system("intersectBed -a " + blastnRes2 + " -b " + altBed + " -wo |awk '{if ($6>90) print}' |cut -f 10 |sort -k1,1 |uniq -d -c |sed 's/ \+/\t/g' > " + outdir + "/query.repeated.hit.txt")
    os.system("intersectBed -a " + blastnRes2 + " -b " + altBed + " -wo |awk '{if ($6>90) print}' |cut -f 4 |sort -k1,1 |uniq -d -c  |sed 's/ \+/\t/g' > " + outdir + "/Col.repeated.hit.txt")
    
    fi = open (outdir + "/Col.repeated.hit.txt","r")
    repHit1 = {}
    #id = ""
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[1].split(".")[0]
        #id = id.split("=")[1]
        repHit1[id] = int(t[0])
        #print id
    fi.close()
    print "\n###\nrepeated hit for Col [mis-split] : ", len(repHit1)
    
    fi = open (outdir + "/query.repeated.hit.txt","r")
    repHit2 = {}
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[1].split(";")[0]
        id = id.split("=")[1]
        repHit2[id] = int(t[0])
    fi.close()
    print "repeated hit for acc [mis-merge] : ", len(repHit2)
    
    fi = open( outdir + "/Col.gene.blastn.intersect.Qry.gene.txt", "r")
    fo1 = open(outdir + "/potential.Qry.missing.gene.txt", "w" )
    fo2 = open(outdir + "/potential.Qry.mis-exon-intron.gene.txt", "w")
    fo3 = open(outdir + "/potential.Qry.mis-split.gene.by.blastn.txt", "w")
    fo4 = open(outdir + "/potential.Qry.mis-merged.gene.by.blastn.txt", "w")
    fo5 = open(outdir + "/potential.Qry.m-vs-m.toBeChecked.by.blastn.txt", "w")
    fo6 = open(outdir + "/futher.check.list","w")
    
    #repHit
    repHitM1 = {} ##col-repeated  mis-split
    repHitM2 = {} ##query-repeated  mis-merge
    k1 = 0
    kM = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        id = t[3].split(".")[0]
        if t[6] == "." :
            fo1.write(line.strip())
            if ungrRef.has_key(id) :
                missingRef[id] = t[0:6]         
                fo1.write("\tungrouped\n")
            else :
                fo1.write("\tgrouped\n")
            kM += 1
        else :            
            refId = t[3].split(".")[0]
            altId = t[9].split(";")[0]
            altId = altId.split("=")[1]
                            
            if float(t[4]) < 90 or float(t[5]) < 90: continue
                
            if repHit1.has_key(refId) and repHit2.has_key(altId):
                
                if repHitM1.has_key(refId) :
                    repHitM1[refId].append(t)
                else :
                    repHitM1[refId] = []
                    repHitM1[refId].append(t)
                gr1 = "ref-grouped"
                gr2 = "qry-grouped"
                if ungrRef.has_key(refId) : 
                    gr1= "ref-ungrouped"
                if ungrAlt.has_key(altId) : 
                    gr2= "alt-ungrouped"
                
                fo5.write("{}\t{}\t{}\n".format(line.strip(), gr1, gr2))  
                #if repHitM2.has_key(altId) :
                #    repHitM2[altId].append(t)
                #else :
                #    repHitM2[altId] = []
                #    repHitM2[altId].append(t)
                    
            elif repHit1.has_key(refId) :
                if repHitM1.has_key(refId) :
                    repHitM1[refId].append(t)
                else :
                    repHitM1[refId] = []
                    repHitM1[refId].append(t)
                                
            elif repHit2.has_key(altId) :
                if repHitM2.has_key(altId) :
                    repHitM2[altId].append(t)
                else :
                    repHitM2[altId] = []
                    repHitM2[altId].append(t)                            
            else :
                if ort.has_key(refId) and (ort[refId].has_key(altId)) :
                    k1 += 1
                    continue
                if ort.has_key(refId) and (not ort[refId].has_key(altId)) :            
                    
                    flag1 = "ref-grouped"
                    flag2 = "alt-grouped"
                    if ungrRef.has_key(refId) : 
                        flag1 = "ref-ungrouped"
                        
                    if ungrAlt.has_key(altId) : 
                        flag2 = "alt-ungrouped"
                        
                    fo2.write("\t".join(t[0:11]))
                                        
                    if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-ortholog"))
                    else :                        
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                        misAnnRef[refId] = 1
                        misAnnAlt[altId] = 1
                elif not ort.has_key(refId) :            
                    #fo2.write(line.strip())
                    flag1 = "ref-grouped"
                    flag2 = "alt-grouped"
                    if ungrRef.has_key(refId) : 
                        flag1 = "ref-ungrouped"
                        
                    if ungrAlt.has_key(altId) : 
                        flag2 = "alt-ungrouped"
                             
                    fo2.write("\t".join(t[0:11]))       
                    
                    if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-ortholog"))
                    else :                        
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                        misAnnRef[refId] = 1
                        misAnnAlt[altId] = 1
                        
                elif float(t[4]) == 100.0 and float(t[5]) == 100.0 and (t[1] != t[7] or t[2] != t[8] ):
                    flag1 = "ref-grouped"
                    flag2 = "alt-grouped"
                    if ungrRef.has_key(refId) : 
                        flag1 = "ref-ungrouped"
                        
                    if ungrAlt.has_key(altId) : 
                        flag2 = "alt-ungrouped"
                        
                    fo2.write("\t".join(t[0:11]))
                    
                    if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-ortholog"))
                    else :                        
                        fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                        misAnnRef[refId] = 1
                        misAnnAlt[altId] = 1
    fi.close()
    print "\n##\npotential missing genes in accession but annotated in Col, ", kM
    print "potential missing genes in accession but annotated in Col, resulted ungrouped Col, ", len(missingRef)
   
    
    for k in sorted(repHitM1.keys()) :       
        if len(repHitM1[k]) == 2 :
            altId1 = repHitM1[k][0][9].split(";")[0]
            altId1 = altId1.split("=")[1]
            altId2 = repHitM1[k][1][9].split(";")[0]
            altId2 = altId2.split("=")[1]
            if ( float(repHitM1[k][0][8]) - float(repHitM1[k][0][7]) ) * 0.9 < float(repHitM1[k][0][-1]) \
                and ( float(repHitM1[k][1][8]) - float(repHitM1[k][1][7]) ) * 0.9 < float(repHitM1[k][1][-1]) and abs(altIdx[altId1] - altIdx[altId2]) == 1 :
                t1 = repHitM1[k][0]
                t2 = repHitM1[k][1]
                #fo3.write("\t".join(t1[0:11]))
                fo3.write("\t".join(t1[0:10]))
                fo3.write("\n")
                fo3.write("\t".join(t2[0:10]))
                fo3.write("\n")
            else :
                fo6.write("\t".join(repHitM1[k][0]))
                fo6.write("\t".join(repHitM1[k][1]))
                fo6.write("\n")
        elif len(repHitM1[k]) > 2:
            flag = 0
        
            for i in range(len(repHitM1[k])) :
                if ( float(repHitM1[k][i][8]) - float(repHitM1[k][i][7]) ) * 0.9 < float(repHitM1[k][0][-1]) :
                    flag = 1
                else :
                    flag = 0
                    break
            if flag == 1 :
                for i in range(len(repHitM1[k])) :        
                    t = repHitM1[k][i]
                    fo3.write("\t".join(t[0:10]))               
                    fo3.write("\n")
            else :
                for i in range(len(repHitM1[k])) :        
                    t = repHitM1[k][i]
                    fo6.write("\t".join(t[0:10]))               
                    fo6.write("\n")
        else :
            fo6.write("\t".join(repHitM1[k][0]))
            fo6.write("\n")
    
    
             
           
    for k in sorted(repHitM2.keys()) :        
        if len(repHitM2[k]) == 2 :
            if ( float(repHitM2[k][0][2]) - float(repHitM2[k][0][1]) ) * 0.9 < float(repHitM2[k][0][-1]) \
                and ( float(repHitM2[k][0][8]) - float(repHitM2[k][0][7]) ) * 0.9 < float(repHitM2[k][0][-1]) \
                and float(repHitM2[k][1][8]) - float(repHitM2[k][1][7]) *0.9 > float(repHitM2[k][1][-1]) :
                ## two genes annotated in one region, may be +/- or same strand (should not be wrong annotation??) 
                t = repHitM2[k][1]
                missingRef[id] = t[0:6]  
                fo1.write("\t".join(t[0:10]))
                fo1.write("\n")
            elif ( float(repHitM2[k][1][2]) - float(repHitM2[k][1][1]) ) * 0.9 < float(repHitM2[k][1][-1]) \
                and ( float(repHitM2[k][1][8]) - float(repHitM2[k][1][7]) ) * 0.9 < float(repHitM2[k][1][-1]) \
                and float(repHitM2[k][0][8]) - float(repHitM2[k][0][7]) *0.9 > float(repHitM2[k][0][-1]) :
                ## two genes annotated in one region, may be +/- or same strand (should not be wrong annotation??)
                t = repHitM2[k][0]
                missingRef[id] = t[0:6]  
                fo1.write("\t".join(t[0:10]))
                fo1.write("\n")
            else  :
                refId1 = repHitM2[k][0][3].split(".")[0]
                refId2 = repHitM2[k][1][3].split(".")[0]
                if abs(refIdx[refId1] - refIdx[refId2]) == 1 and  repHitM2[k][0][1] != repHitM2[k][1][1]:
                    t1 = repHitM2[k][0]
                    t2 = repHitM2[k][1]
                    #if k == "ATAN1-1G16270" :
                    #    print t1,"\n"
                    #    print t2,"\n"
                    fo4.write("\t".join(t1[0:10]))
                    fo4.write("\n")
                    fo4.write("\t".join(t2[0:10]))
                    fo4.write("\n")
                else :
                    if  float(repHitM2[k][0][4]) == 100.0 and float(repHitM2[k][0][5]) == 100.0 and ( repHitM2[k][0][1] != repHitM2[k][0][7] or repHitM2[k][0][2] != repHitM2[k][0][8] ) :
                        flag1 = "ref-grouped"
                        flag2 = "alt-grouped"
                        flag3 = "not-orthlog"                        
                        refId = repHitM2[k][0][3].split(".")[0]
                        altId = k
                        t = repHitM2[k][0]
                        fo2.write("\t".join(t[0:10]))
                        if ort.has_key(refId) and ort[refId].has_key(k) :
                            flag3 = "orthlog"
                        if ungrRef.has_key(refId) : 
                            flag1 = "ref-ungrouped"
                            
                        if ungrAlt.has_key(altId) : 
                            flag2 = "alt-ungrouped"
                               
                        if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                            fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-orthloog"))
                        else :                        
                            fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                            misAnnRef[refId] = 1
                            misAnnAlt[altId] = 1
                        #fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,flag3))       
                    elif  float(repHitM2[k][1][4]) == 100.0 and float(repHitM2[k][1][5]) == 100.0 and ( repHitM2[k][1][1] != repHitM2[k][1][7] or repHitM2[k][1][2] != repHitM2[k][1][8] ) :
                        flag1 = "ref-grouped"
                        flag2 = "alt-grouped"
                        flag3 = "not-orthlog"                        
                        refId = repHitM2[k][0][3].split(".")[0]
                        altId = k
                        t = repHitM2[k][1]
                        fo2.write("\t".join(t[0:10]))
                        if ort.has_key(refId) and ort[refId].has_key(k) :
                            flag3 = "orthlog"
                        if ungrRef.has_key(refId) : 
                            flag1 = "ref-ungrouped"
                            
                        if ungrAlt.has_key(altId) : 
                            flag2 = "alt-ungrouped"
                            
                        if blastp.has_key(altId + refId) and blastp[altId + refId][0] > 80 and blastp[altId + refId][1] > 0.80 and blastp[altId + refId][2] > 0.80:
                            fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"missing-ortholog"))
                        else :                        
                            fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,"not-ortholog"))
                            misAnnRef[refId] = 1
                            misAnnAlt[altId] = 1   
                        #fo2.write("\t{}\t{}\t{}\n".format(flag1,flag2,flag3))           
                                                
                    
        elif len(repHitM2[k]) > 2:
            flag = 0
            for i in range(len(repHitM2[k])) : 
                refId = repHitM2[k][i][3].split(".")[0]
                #print refId
                if ort.has_key(refId) and ort[refId].has_key(k) and (float(repHitM2[k][i][8]) - float(repHitM2[k][i][7]) ) * 0.8 < float(repHitM2[k][i][-1]):
                    flag = 1  ### multiple vs one, some deleleted in query genome
            if int(repHitM2[k][1][1]) < int(repHitM2[k][0][2]) and  int(repHitM2[k][1][1]) - int(repHitM2[k][0][1]) <  int(repHitM2[k][0][2]) - int(repHitM2[k][1][1]) :
                flag = 1
            
            if flag == 0 :
                for i in range(len(repHitM2[k])) :
                    t = repHitM2[k][i]
                    fo4.write("\t".join(t[0:10]))
                    fo4.write("\n")
            else :
                for i in range(len(repHitM2[k])) :
                    t = repHitM2[k][i]
                    fo6.write("\t".join(t[0:10]))
                    fo6.write("\n")
        else :
            fo6.write("\t".join(repHitM2[k][0]))
            fo6.write("\n")
    fo1.close()
    fo2.close()
    fo3.close()
    fo4.close()
    fo5.close()
    fo6.close()
    print "\n#false mis-ex as they are orthlog: ",  k1
    print "\n"
    return [missingRef,misAnnRef,misAnnAlt,missAssRef]
            
def getAltSpeGene(ungrAlt,altPos,outfile):  
    fo = open(outfile,"w")
    #print len(ungrAlt.keys())
    for k in sorted(ungrAlt.keys()):
        [chrom,start,end] = altPos[k]
        fo.write("{}\t{}\t{}\t{}\n".format(chrom,start,end,k))
        #sys.exit()
    fo.close()
    

def getProtLeng2(protFas):
    fx = open(protFas,"r")
    leng = 0
    l = {}
    id = ""
    while True :
        line = fx.readline()
        if not line : break
        if line[0]==">" :
            if id != "" :
                id = id.split("|")[1]
                if re.search("evm",id) :
                    id = re.sub("model","TU",id)
                else :
                    id = id.split(".")[0]                
                l[id] = leng
                leng = 0            
            t = line.strip().split()
            id = t[0]
        else :
            leng += len(line.strip())
    id = id.split("|")[1]
    if re.search("evm",id) :
        id = re.sub("model","TU",id)
    else :
        id = id.split(".")[0]            
    l[id] = leng
    fx.close()
    return l

def getProtLeng(protFas, flag):
    l = {}
    if not os.path.isfile(protFas) :
        print "no such file ", protFas
        sys.exit()
    for seq_record in SeqIO.parse(protFas, "fasta"):
        #print(seq_record.id)
        #print(repr(seq_record.seq))
        #print(len(seq_record))
        if flag == 1 :
            l[seq_record.id.split(".")[0]] = len(seq_record)
        else :
            l[seq_record.id] = len(seq_record)
    return l




def getSpecGene(groupFile,protRef,protAlt):
    fx = open(groupFile,"r")
    gr = {}
    ort = {}
    ref = "AT"
    while True :
        line = fx.readline()
        if not line : break       
        if not re.search(ref, line) or not re.search("evm",line) : continue 
        t = line.strip().split()
        refIds = []
        altIds = []        
        for k in range(len(t)) :
            if k == 0 : continue
            if re.search(ref,t[k]) :
                id = t[k].split(".")[0]
                refIds.append(id)
                gr[id]=1
            else :            
                #id = t[k].split(".")[0]                
                id = t[k]
                altIds.append(id)            
                gr[id]=1
        for k in refIds :
            ort[k]={}
            for j in altIds :
                ort[k][j]=1                        
    fx.close()    
    ungrRef = {}
    ungrAlt = {}
    
    for k in sorted(protRef.keys()) :
        if not gr.has_key(k) :
            ungrRef[k]=1
    for k in sorted(protAlt.keys()) :
        if not gr.has_key(k) :
            ungrAlt[k]=1            
    print "ungrouped ref gene ", len(ungrRef)
    print "ungrouped alt gene ", len(ungrAlt)
    return [ungrRef,ungrAlt,ort]
        
    
if __name__ == "__main__":
   main(sys.argv[1:])    
   
   
   