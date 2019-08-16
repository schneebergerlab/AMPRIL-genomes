#!/usr/bin/env python
# encoding: utf-8

'''
Created on Aug 02, 2018
build pangenome based on protein-coding genes orthologs resulting from orthoFinder
@author: jiao
'''
import sys
import os
import argparse
import datetime
from scipy.special import comb, perm
from itertools import combinations, permutations


def main(argv):
    parser = argparse.ArgumentParser(description="")

    parser.add_argument('-g','--group', required=True, help="gene family output file from orthofinder")
    #parser.add_argument('-f','--faDir', required=True, help="all fasta file dir")    
    parser.add_argument('-o','--outdir', default=False ,help="output directory, default: current working directory ")
    
    args = parser.parse_args()
    
    startTime = datetime.datetime.now()
    print "Starting ======  ",startTime.strftime('%Y-%m-%d %H:%M:%S'), "  ======"
    #faDir = args.faDir
    outdir = args.outdir
    grFile = args.group
    
    

    
    accs = ["An-1","C24","Col","Cvi","Eri","Kyo","Ler","Sha"]
    
    startTime = datetime.datetime.now()
    print "Starting ======  ",startTime.strftime('%Y-%m-%d %H:%M:%S'), "  ======"
    if not os.path.isdir(outdir) :
        os.system("mkdir -p " + outdir)
        
   
    arr = []

    for i in range(len(accs)) :
        arr.append(i)
       
    #st = [[],[]] #core consensus
    [n, genes,groups] = getGens(grFile,accs)
    
    results = []
    for i in range(len(accs)) :
        num = i + 1
        result = runJob(n, genes, groups,accs, arr, num)
        results.append(result)    
    
    outFile1 = outdir + "/pan-genome.gene.core.stats"
    outFile2 = outdir + "/pan-genome.gene.conensus.stats"
    outFile3 = outdir + "/pan-genome.newGenes.stats"        
    fo1 = open(outFile1, "w")
    fo2 = open(outFile2, "w")
    fo3 = open(outFile3, "w")        
    for i in range(len(accs)) :
        num = i + 1
        print "# number of genomes: ", str(num)
        coms = list(combinations(arr,num))
        print "    number of combinations: ", len(coms), " x len(accs)"
           
        stat1 = results[i][0]
        stat2 = results[i][1]
        fo1.write("{}\t{}\n".format(i, "\t".join(str(x) for x in stat1)))
        fo2.write("{}\t{}\n".format(i, "\t".join(str(x) for x in stat2)))
        if i > 0 :
            stat3 = results[i][2]
            fo3.write("{}\t{}\n".format(i, "\t".join(str(x) for x in stat3)))    
                
    fo1.close()
    fo2.close()
    fo3.close()

    endTime = datetime.datetime.now()
    print "Finished ... Time elapsed: ",(endTime - startTime).seconds, " s"
    
def getGens (inFile,accs): 
    genes = {}
    groups = []
    for acc in accs :
        genes[acc] = []
    n = 0
    fi = open(inFile,"r")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        for i in range(len(accs)) :
            genes[accs[i]].append(t[i+1])
        groups.append(t[0])
        n += 1
        #print genes["An-1"][3333]
        #sys.exit()
    fi.close()
    #print len(genes[accs[3]])
    return [n, genes, groups]
    
def runJob(n, genes, groups, accs, arr, num):
 
    stat1 = []
    stat2 = []
    stat3 = []
    coms = list(combinations(arr,num))            
    for k in coms : #[0,1,2],[1,2,3] [2,3,4]
        comAccs = []
        for j in k :
            comAccs.append(accs[j])  
                        
        core = getCore(genes, n, comAccs)
        stat1 = stat1 + core
        print "\t\t combination: ", k , comAccs, "\t core :" , core
                   
        cons = getCons(genes, groups, n, comAccs)
        stat2 = stat2 + cons
        print "\t\t combination: ", k , comAccs, "\t consenus :" , cons
        
        newG = getNewGenes(genes, groups,n, comAccs)
        stat3 = stat3 + newG
    print [stat1, stat2, stat3]
    return [stat1,stat2, stat3]

def getCons(genes, groups, n, subAccs) :
    cons = []
    #print "for subAccs", subAccs
    for acc1 in subAccs :   
        geneNum = 0     
        useGr = {}  
        for i in range(n) :
            if genes[acc1][i] != "-----" :
                geneNum = geneNum + len(genes[acc1][i].split(","))
                if not useGr.has_key(groups[i]) :
                    useGr[groups[i]] = 1
                  
            #maxN = 0        
        for acc2 in subAccs :
            if acc2 == acc1 : continue
            for i in range(n) :
                if genes[acc1][i] != "-----"  : continue                             
                #if genes[acc2][i] != "-----" and len(genes[acc2][i].split(",")) > maxN :
                if genes[acc2][i] != "-----" :
                    #maxN = len(genes[acc2][i].split(","))                                        
                    if not useGr.has_key(groups[i]) :
                        geneNum += len(genes[acc2][i].split(","))
                        useGr[groups[i]] = 1
        print subAccs, geneNum
        cons.append(geneNum)            
    #avg = int(sum(cons)/len(cons))
    #return avg
    return cons

def getCore(genes, n , subAccs):
    core = []
    for acc1 in subAccs :   
        geneNum = 0     
        for i in range(n) :
            flag = 1
            for acc2 in subAccs :
                if genes[acc2][i] == "-----"  : flag = 0                                
            if flag == 1 :
                geneNum = geneNum + len(genes[acc1][i].split(","))
                
        core.append(geneNum)            
    #avg = int(sum(core)/len(core))
    #return avg
    return core

def getNewGenes(genes, groups,n, subAccs):
    newGenes = []
    #print "for subAccs", subAccs
    for acc1 in subAccs :   
        newG = 0        
        for i in range(n) :
            if genes[acc1][i] != "-----" :
                newG += 1
        useGr = {}
        for acc2 in subAccs :
            if acc2 == acc1 : continue
            for i in range(n) :
                if genes[acc1][i] == "-----"  : continue                             
                #if genes[acc2][i] != "-----" and len(genes[acc2][i].split(",")) > maxN :
                if genes[acc2][i] != "-----" :
                    #maxN = len(genes[acc2][i].split(","))  
                    if not useGr.has_key(groups[i]) :                                                                             
                        newG -= 1
                        useGr[groups[i]] = 1                     
        #print subAccs, newG
        newGenes.append(newG)            
    #avg = int(sum(cons)/len(cons))
    #return avg
    return newGenes
    

if __name__ == "__main__":
   main(sys.argv[1:])

