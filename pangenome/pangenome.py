#!/usr/bin/env python
# encoding: utf-8

import sys
import os
import datetime
import getopt
from scipy.special import comb, perm
from itertools import combinations, permutations
#list(combinations([1, 2, 3], 2))
#[(1, 2), (1, 3), (2, 3)]

def main(argv):
    wgadir = ""
    outdir = ""
    gBedDir = ""
    try:
        opts, args = getopt.getopt(argv,"w:o:g:",["wgadir=","outdir=","gBed=",])
    except getopt.GetoptError:
        print 'wga.pangenome.py -w <wgadir> -o <outdir>  -g <gBed> '
        sys.exit(2)
    if len(opts) == 0 :
        print 'wga.pangenome.py -w <wgadir> -o <outdir>  -g <gBed> '
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print 'wga.pangenome.py -w <wgadir> -o <outdir>  -g <gBed> '
            sys.exit()
        elif opt in ("-w", "--wgadir"):
            wgadir = arg
        elif opt in ("-o", "--outdir"):
            outdir = arg
        elif opt in ("-g", "--gBed") :
            gBedDir = arg 
    
    startTime = datetime.datetime.now()
    print "Starting ======  ",startTime.strftime('%Y-%m-%d %H:%M:%S'), "  ======"
    if not os.path.isdir(outdir + "/tmp") :
        os.system("mkdir -p " + outdir + "/tmp")
    accs = ["An-1","C24","Col","Cvi","Eri","Kyo","Ler","Sha"]
    arr = []
    chrBeds = {}
    for i in range(len(accs)) :
        arr.append(i)
        genomeFile = gBedDir + "/" + accs[i] + ".leng.txt"
        chrBeds[accs[i]]= genomeFile
    
    #st = [[],[]] #core consensus/pan
    results = []
    for i in range(len(accs)) :
        num = i + 1
        result = runJob(arr, num, chrBeds, accs, wgadir, outdir)
        results.append(result)
    
    
    outFile1 = outdir + "/pan-genome.wga.core.stats"
    outFile2 = outdir + "/pan-genome.wga.conensus.stats"
    outFile3 = outdir + "/pan-genome.wga.newseq.stats"     
    fo1 = open(outFile1, "w")
    fo2 = open(outFile2, "w")
    fo3 = open(outFile3, "w")          
    for i in range(len(accs)) :
        num = i + 1
        print "# number of genomes: ", str(num)
        coms = list(combinations(arr,num))
        print " the number N of independent measurements: ", len(coms)*len(coms[0])
        #stat1 = results[i].get()[0]
        #stat2 = results[i].get()[1]        
        stat1 = results[i][0]
        stat2 = results[i][1]
        stat3 = results[i][2]
        fo1.write("{}\t{}\n".format(i, "\t".join(str(x) for x in stat1)))
        fo2.write("{}\t{}\n".format(i, "\t".join(str(x) for x in stat2)))        
        fo3.write("{}\t{}\n".format(i, "\t".join(str(x) for x in stat3)))
    fo1.close()
    fo2.close()
    fo3.close()

    endTime = datetime.datetime.now()
    print "Finished ... Time elapsed: ",(endTime - startTime).seconds, " s"
    
def runJob(arr, num, chrBeds, accs, wgadir, outdir):
 
    stat1 = []
    stat2 = []
    stat3 = []
    coms = list(combinations(arr,num))
    if num == 1 :
        stat1 = []
        for j in range(len(accs)) :
            t = getLen(chrBeds[accs[j]])
            stat1.append(t)
        stat2 = stat1
        stat3 = stat1
    else :            
        for k in coms : #[0,1,2],[1,2,3] [2,3,4]
            comAccs = []
            for j in k :
                comAccs.append(accs[j])  
            
            st = getCore(comAccs, wgadir, chrBeds, outdir)
            stat1 = stat1 + st
            print "\t\t combinations: ", k , comAccs, "\t core :" , st
            
            st = getNewSeq(comAccs, wgadir, chrBeds, outdir)
            stat3 = stat3 + st
            print "\t\t combinations: ", k , comAccs, "\t new seq :" , st            
            
            st = []
            for j in range(len(comAccs)) :                           
                x = getCons(comAccs[j],comAccs, wgadir, chrBeds, outdir, 0 )
                st.append(x)                
            stat2 = stat2 + st
            print "\t\t combination: ", k , comAccs, "\t consenus :" , st
            
    print [stat1, stat2, stat3]
    return [stat1,stat2, stat3]

def getCore(subAccs, wgadir, chrBeds, outdir) :
    core = []
    for acc1 in subAccs :
        delBedFiles = []
        for acc2 in subAccs :
            if acc1== acc2 : continue            
            delBedFile = os.path.join(wgadir, acc1, acc2, acc2 + ".del.bed")
            delBedFiles.append(delBedFile)
        outFile = outdir + "/tmp/" + acc1 + ".com-" + str(len(subAccs)) + ".core.bed"            
        cmd1 = "cat " + " ".join(delBedFiles) + " | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i -  | "
        cmd2 = "bedtools complement -i - -g " + chrBeds[acc1] + " > " + outFile
        print cmd1 + cmd2  
        os.system(cmd1 + cmd2)
        t = getLen2(outFile)
        core.append(t)
    #avg = int(sum(core)/len(core))
    #return avg
    return core

def getCons2(subAccs, wgadir, chrBeds, outdir):
    cons = []
    for acc1 in subAccs :
        con = getLen(chrBeds[acc1])        
        for acc2 in subAccs :
            if acc1== acc2 : continue            
            insBedFile = os.path.join(wgadir, acc1, acc2, acc2 + ".ins.bed")
            for acc3 in subAccs :
                if acc3 == acc2 : continue
                wgaFile = os.path.join(wgadir, acc2, acc3, acc3 + ".wga.block.txt")
                outFile = outdir + "/tmp/tmp." + acc1 + "." + acc2 + "." + acc3 + ".bed"
                cmd = "bedtools subtract -a " + insBedFile + " -b " + wgaFile + " > " + outFile  
                print cmd               
                os.system(cmd)
                insBedFile = outFile
            t = getLen2(insBedFile)
            con = con + t        
        cons.append(con)
    avg = int(sum(cons)/len(cons))
    return avg
    return

def getCons(acc1, subAccs, wgadir,chrBeds, outdir, x):    
    if len(subAccs) == 1 :
        acc = subAccs[0]
        x = x + getLen(chrBeds[acc])
        #y = getLen(chrBeds[acc])
        print "con leng" , x
        #print "new seq", y
        return x
    else :
        inFile = chrBeds[acc1]
        genomeBed = outdir + "/tmp/" + acc1 + ".genome.bed" 
        getGenomeBed(inFile,genomeBed)
        insBedFile = genomeBed
        for acc2 in subAccs :
            if acc2 == acc1  : continue  
            wgaFile = os.path.join(wgadir, acc1, acc2, acc2 + ".wga.block.txt")
            outFile = outdir + "/tmp/tmp." + acc1 + "." + acc2  + ".bed"
            cmd = "bedtools subtract -a " + insBedFile + " -b " + wgaFile + " > " + outFile  
            print cmd               
            os.system(cmd)
            insBedFile = outFile
        #y = getLen2(insBedFile)
        x = x + getLen2(insBedFile)        
        if len(subAccs) >= 2 :
            return getCons(subAccs[1], subAccs[1:], wgadir, chrBeds, outdir, x)        
        
def getNewSeq(subAccs, wgadir, chrBeds, outdir): 
    newSeq = []
    for acc1 in subAccs :
        wgaBedFiles = []
        for acc2 in subAccs :
            if acc1== acc2 : continue            
            wgaBedFile = os.path.join(wgadir, acc1, acc2, acc2 + ".wga.block.txt")
            wgaBedFiles.append(wgaBedFile)
        outFile = outdir + "/tmp/" + acc1 + ".com-" + str(len(subAccs)) + ".wga.bed"            
        cmd1 = "cat " + " ".join(wgaBedFiles) + " | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i -  | "
        cmd2 = "bedtools complement -i - -g " + chrBeds[acc1] + " > " + outFile
        print cmd1 + cmd2  
        os.system(cmd1 + cmd2)
        t = getLen2(outFile)
        newSeq.append(t)
    #avg = int(sum(core)/len(core))
    #return avg
    return newSeq

def getGenomeBed(inFile,outFile):
    fi = open(inFile,"r")
    fo = open(outFile,"w")
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        fo.write("{}\t{}\t{}\n".format(t[0],0,int(t[1])-1))
    fi.close()
    fo.close()

def getLen(inFile) : #genome len file
    fi = open(inFile,"r")
    total = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        total+=int(t[1])
    return total

def getLen2(inFile) : #bed
    fi = open(inFile,"r")
    total = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split("\t")
        total+=int(t[2]) - int(t[1])
    return total
if __name__ == "__main__":
   main(sys.argv[1:])






 
                