#!/usr/bin/env python
# encoding: utf-8
'''
Created on Jul 26, 2017

@author: jiao
'''
import glob
import os
import re
import time
import sys
import subprocess

def splitProt (fastaDir,chunkSize,outdir):
    fastaFiles = glob.glob(fastaDir + "/*.fa*")
    n = 0
    chunk = 0
    file = "prot.chunk." + str(chunk) + ".fa" 
    fo = open(os.path.join(outdir,file),"w")
    #print "chunkSize",chunkSize
    for fasta in fastaFiles :                
        fi = open(fasta,"r")        
        while True:
            line = fi.readline()
            if not line : break
            if line[0] == ">" :
                #print line
                if n == int(chunkSize) :
                    fo.close()
                    chunk += 1
                    file = "prot.chunk." + str(chunk) + ".fa" 
                    fo = open(os.path.join(outdir,file),"w")
                    fo.write(line)
                    n = 1                    
                else :
                    n += 1
                    #print n
                    fo.write(line)
            else :
                fo.write(line)
        fi.close()
    fo.close()
        

def parserConfig (cfgFile):    
    cfg = {}
    fi = open(cfgFile,"r")
    while True :
        line = fi.readline()
        if not line : break
        line = line.strip()
        if re.search(r"#",line) : continue
        if not re.search(r'=',line): continue
        t = line.split("=")
        t[0] = t[0].strip()
        t[1] = t[1].strip()        
        cfg[t[0]] = t[1]
        #print t[0],t[1]
        if re.search(r"Outdir",t[0]) :
            if not os.path.isdir(t[1]) :
                os.system("mkdir -p " + t[1])
                print "mkdir ",t[0]
        if t[0]=="protAligner" :
            outdir = os.path.join(cfg["protOutdir"],t[1])
            if not os.path.isdir(outdir) :
                os.system("mkdir -p " + outdir)  
    return cfg

def runHisatStringtie (name,read1,read2,ref,outdir,thread,out):    
    sumfile = os.path.join(outdir, name + ".summary")
    samout = os.path.join(outdir, name + ".sam")
    bamout = os.path.join(outdir, name + ".bam")
    srtout = os.path.join(outdir, name + ".srt.bam")
    gtfout = os.path.join(outdir, name + ".gtf")
    merlist =os.path.join(outdir,  "mergelist.txt")
    #if (os.path.isfile(merlist)) : os.system("rm " + merlist)
    try :
        fo = open(out,"w")
        print "echo file", out
    except :
        print "cannot open file ", out
        sys.exit(2)    
    #os.system("hisat2 --dta --summary-file " + sumfile + " -p " + thread + " -x " + ref + " -1 " + read1 + " -2 " + read2 + " -S " + samout)
    if read2 == "None" :
        fo.write("hisat2 --dta --summary-file " + sumfile + " -p " + str(thread) + " -x " + ref + " -U " + read1 + " -S " + samout)
    else :
        fo.write("hisat2 --dta --summary-file " + sumfile + " -p " + str(thread) + " -x " + ref + " -1 " + read1 + " -2 " + read2 + " -S " + samout)
    fo.write("\n")
    #os.system("samtools view -bS " + samout + " -o " + bamout)
    fo.write("samtools view -bS " + samout + " -o " + bamout)
    fo.write("\n")
    #os.system("samtools sort -@ 10 " + " -o " + srtout + " " + bamout)
    fo.write("samtools sort -@ 10 " + " -o " + srtout + " " + bamout)
    fo.write("\n")
    #stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR204916_chrX.gtf –l ERR204916 ERR204916_chrX.bam
    #os.system("stringtie " + srtout + " -p " + str(thread) + " -o " + gtfout )
    #fo.write("\n")
    #stringtie --merge -p 8 -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf chrX_data/mergelist.txt
    fo.write("stringtie " + srtout + " -p " + str(thread) + " -o " + gtfout)
    fo.write("\n")
    #os.system("echo " + gtfout + " >> " + merlist)    
    #os.system("rm " + samout + " " + bamout)
    fo.write("echo " + gtfout + " >> " + merlist)
    fo.write("\n")
    fo.write("rm " + samout + " " + bamout)
    fo.write("\n")
    fo.close()

def splitGenome (genome,chunkSize,outdir):
    fi = open(genome,"r")
    chunk = 0
    cumSize = 0
    seq = ""
    while True :
        line = fi.readline()
        if not line : break
        line.strip()
        if line[0] == ">" :
            if cumSize > int(chunkSize) :
                chunk += 1
                out = os.path.join(outdir,"genome.chunk." + str(chunk) + ".fa")
                fo = open(out,"w")
                fo.write(seq)
                fo.close()
                seq = line
                cumSize = 0
            else :
                seq = seq + line
        else :
            cumSize += len(line)
            seq = seq + line
            
    if cumSize > 0 :
        chunk += 1
        out = os.path.join(outdir, "genome.chunk." + str(chunk) + ".fa")
        fo = open(out,"w")
        fo.write(seq)
        fo.close()
        
def striGTF2augHint (gtfFile,out):
    fi = open(gtfFile,"r")
    fo = open(out,"w")
    while True :
        line = fi.readline()
        if not line : break
        if line[0]=="#" : continue
        t = line.strip().split("\t")
        tt = t[8].split()
        t[2] = "ep"
        t[1] = "b2h"
        id = tt[3][1:int(len(tt[3]))-2]
        #print id     
        #sys.exit()            
        line = "\t".join(t[0:7])
        line = line + "\tgrp=" + id + ";pri=4\;src=W\n"
        fo.write(line)
    fo.close()
    fi.close()



def jobWait(jobs):
    jobDone = {}    
    print "len(jobs)", len(jobs)
    while True :        
        for job in jobs :
            if not os.path.isfile(job) : 
                #time.sleep(5)
                continue
            #time.sleep(1)                 
            fi = open(job,"r")
            while True :
                line = fi.readline()
                if not line : break
                if re.search("Successfully", line) :
                    jobDone[job] = 1
                    break                                
        if len(jobDone.keys()) == len(jobs) :
            break
        else :
            time.sleep(5)
    print "all job finished\n"
    return

def geneStats (gff3,out):
    
    return
    
    
    