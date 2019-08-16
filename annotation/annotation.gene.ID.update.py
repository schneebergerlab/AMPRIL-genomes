#!/usr/bin/env python
# encoding: utf-8


import re
import sys
import os
import glob
import getopt


def main(argv):
    gffI = ""
    gffO = ""
    outdir = ""    
    acc = ""
    ver = ""
    genome = ""
    try:
        opts, args = getopt.getopt(argv,"i:o:n:a:v:g:",["gffI=","gffO=","gffN=","acc=","ver=","ref="]) 
    except getopt.GetoptError:
        print 'annotation.gene.ID.update.py -i <gffI> -o <outdir> -n <gffN> -a <acc> -v <ver> -g <ref> '
        sys.exit(2)
    if len(opts) == 0 :
        print 'annotation.gene.ID.update.py -i <gffI> -o <outdir> -n <gffN> -a <acc> -v <ver> -g <ref> '
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print 'annotation.gene.ID.update.py -i <gffI> -o <outdir> -n <gffN> -a <acc> -v <ver> -g <ref> '
            sys.exit()
        elif opt in ("-i", "--gffI"):
            gffI = arg           
        elif opt in ("-o", "--gffO"):
            outdir = arg
        elif opt in ("-n", "--gffN"):
            gffN = arg
        elif opt in ("-a", "--acc"):
            acc = arg
        elif opt in ("-v", "--ver"):
            ver = arg
        elif opt in ("-g", "--ref"):
            genome = arg
            
    pre = acc.upper()
    if acc=="An-1" : pre = "AN"
    
    fi = open(gffI,"r")
    genes1 = {}  ##chr1-5
    genesStart1 = {}
    genes2 = {} ## unachored
    genesStart2 = {}
    id = ""     
    m = 0    
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split()                            
        if t[2] == "gene" :
            m+=1
            id = t[8].split(";")[0]
            id = id.split("=")[1]
            if re.search("chr", line) :                
                if not genes1.has_key(t[0]) :
                    genes1[t[0]] = {}
                    genes1[t[0]][id] = line
                    genesStart1[t[0]]={}
                    genesStart1[t[0]][id] = t[3]
                else :                
                    genes1[t[0]][id] = line
                    genesStart1[t[0]][id] = t[3]
            else :
                if not genes2.has_key(t[0]) :
                    genes2[t[0]] = {}
                    genes2[t[0]][id] = line
                    genesStart2[t[0]] = {}
                    genesStart2[t[0]][id] = t[3]
                else :                
                    genes2[t[0]][id] = line
                    genesStart2[t[0]][id] = t[3]
        else :                    
            if re.search("chr",line) :
                #if not genes1[t[0]].has_key(id) : 
                #    print line
                #    sys.exit()                    
                genes1[t[0]][id] = genes1[t[0]][id] + line
            else :
                genes2[t[0]][id] = genes2[t[0]][id] + line
            
    fi.close()
    print "protein-coding genes:\t", m
     
    fi = open(gffN,"r")    
    flag = 0    
    m1 = 0
    m2 = 0
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split()                            
        if re.search("Infer",t[1]) :
            m1+=1
            id = t[8].split(";")[0]
            id = id.split("=")[1]
            if re.search("chr", line) :                
                if not genes1.has_key(t[0]) :
                    genes1[t[0]] = {}
                    genes1[t[0]][id] = line
                    genesStart1[t[0]] = {}
                    genesStart1[t[0]][id] = t[3]
                else :                
                    genes1[t[0]][id] = line
                    genesStart1[t[0]][id] = t[3]
            else :
                if not genes2.has_key(t[0]) :
                    genes2[t[0]] = {}
                    genes2[t[0]][id] = line
                    genesStart2[t[0]] = {}
                    genesStart2[t[0]][id] = t[3]
                else :                
                    genes2[t[0]][id] = line
                    genesStart2[t[0]][id] = t[3]
                
        elif t[2] == "gene"  : 
            if not re.search("protein_coding",t[8]) :    
                flag = 1  
                m2+=1              
            else :
                flag = 0                
                continue
            id = t[8].split(";")[0]
            id = id.split("=")[1]
            if re.search("chr", line) :                            
                if not genes1.has_key(t[0]) :
                    genes1[t[0]] = {}
                    genes1[t[0]][id] = line
                    genesStart1[t[0]] = {}
                    genesStart1[t[0]][id] = t[3]
                else :                
                    genes1[t[0]][id] = line
                    genesStart1[t[0]][id] = t[3]
            else :                
                if not genes2.has_key(t[0]) :
                    genes2[t[0]] = {}
                    genes2[t[0]][id] = line
                    genesStart2[t[0]] = {}
                    genesStart2[t[0]][id] = t[3]
                else :                
                    genes2[t[0]][id] = line
                    genesStart2[t[0]][id] = t[3]
        else :
            if flag == 0 : continue
            if re.search("chr",line) :                 
                genes1[t[0]][id] = genes1[t[0]][id] + line
            else :
                genes2[t[0]][id] = genes2[t[0]][id] + line
        
        
    fi.close()
    print "non-coding genes:\t", m1
    print "TE-related genes:\t", m2
    
    tmpFile = outdir + "/genes.tmp.gff"
    fo = open(tmpFile,"w")
    for chrom in sorted(genes1.keys()) :
        genes = genes1[chrom]    
        geneStart = genesStart1[chrom]
        for k in sorted(geneStart.items(), key=lambda x: int(x[1])) :
            id = k[0]
            #print k
            #sys.exit()
            
            lines = genes1[chrom][id]
            fo.write(lines)
    
    for chrom in sorted(genes2.keys()) :
        genes = genes2[chrom]    
        geneStart = genesStart2[chrom]
        for k in sorted(geneStart.items(), key=lambda x: int(x[1])) :
            id = k[0]
            lines = genes2[chrom][id]
            fo.write(lines)                
    fo.close()
    
    fi = open(tmpFile,"r")
    outFile1 = outdir + "/" + acc + ".genes.annotation." + ver +".gff"
    outFile2 = outdir + "/" + acc + ".protein-coding.genes."+ ver +".gff"
    outFile3 = outdir + "/" + acc + ".gene.ID.translation." + ver + '.txt'
    fo1 = open(outFile1,"w")
    fo2 = open(outFile2,"w")
    fo3 = open(outFile3,"w")
    gID = ""
    mID = 0
    chrom = ""
    idx = 10010
    idxUn = ""
    
    cdsID = 0
    flag = "prot"
    while True :
        line = fi.readline()
        if not line : break
        t = line.strip().split()
        if not re.search("chr",t[0]) :
            if re.search("Infer",t[1]) :
                if idxUn == "" :
                    idxUn = 10010
                else :
                    idxUn += 10
                infor = "ID=AT" + pre + "-UG" + str(idxUn)
                infor = infor + ";Note=" + t[2]
                outline ="\t".join(t[0:8]) + "\t" + infor + "\n"
                fo1.write(outline) 
                id = t[8].split(";")[0]
                id = id.split("=")[1] 
                newId =  "AT" + pre + "-UG" + str(idxUn)
                fo3.write("{}\t{}\n".format(id,newId))
            else :                
                if t[2] == "gene" :
                    note = t[8].split(";")[-1]
                    note = note.split("=")[1]
                    if re.search("transp",note) : 
                        flag = "TE"
                    else :
                        flag = "prot"
                    if idxUn == "" :                    
                        idxUn = 10010
                    else :
                        idxUn += 10
                    infor = "ID=AT" + pre + "-UG" + str(idxUn)
                    if flag == "prot" :
                        infor = infor + ";Note=protein_coding_gene"
                    else :
                        infor = infor + ";Note=transposable_element_gene"  
                    mID = 0
                    id = t[8].split(";")[0]
                    id = id.split("=")[1] 
                    newId =  "AT" + pre + "-UG" + str(idxUn)
                    fo3.write("{}\t{}\n".format(id,newId))
                    
                elif t[2] == "mRNA" :
                    mID += 1
                    infor = "ID=AT" + pre + "-UG" + str(idxUn) + "." + str(mID) + ";Parent=AT" + pre + "-UG" + str(idxUn)
                    
                elif t[2] =="exon" :
                    nn = t[8].split(";")[0]
                    nn = nn.split("=")[1]
                    tt = nn.split(".")
                    mm = tt[-1]
                    infor = "ID=AT" + pre + "-UG" + str(idxUn) + "." + str(mID) + "." + mm 
                    cdsID = re.sub(r"exon","", mm)
                    infor = infor + ";Parent=AT" + pre + "-UG" + str(idxUn) + "." + str(mID)
                elif t[2] =="CDS" :
                    nn = t[8].split(";")[0]
                    nn = nn.split("=")[1]
                    tt = nn.split(".")
                    #cdsId = tt[-1]
                    #infor = "ID=AT" + pre + "-" + chrom + "G" + str(idxUn) + "." + str(mID) +  ".cds" + str(cdsID) 
                    infor = "ID=AT" + pre + "-UG" + str(idxUn) + "." + str(mID) +  ".cds" + str(cdsID)
                    infor = infor + ";Parent=AT" + pre + "-UG" + str(idxUn) + "." + str(mID)                    
                outline ="\t".join(t[0:8]) + "\t" + infor + "\n"
                fo1.write(outline)  
                if flag == "prot" : 
                    fo2.write(outline)
                        
            continue
       
        if re.search("Infer",t[1]) :
            if not chrom :
                chrom = t[0][3]
                idx = 10010
            elif chrom != t[0][3] :
                chrom = t[0][3]
                idx = 10010
            else :
                idx += 10
            infor = "ID=AT" + pre + "-" + chrom + "G" + str(idx)
            infor = infor + ";Note=" + t[2]
            outline ="\t".join(t[0:8]) + "\t" + infor + "\n"
            fo1.write(outline) 
            id = t[8].split(";")[0]
            id = id.split("=")[1] 
            newId =  "AT" + pre + "-" + chrom + "G" + str(idx)
            fo3.write("{}\t{}\n".format(id,newId))
        else :
            
            if t[2] == "gene" :
                note = t[8].split(";")[-1]
                note = note.split("=")[1]
                if re.search("transp",note) : 
                    flag = "TE"
                else :
                    flag = "prot"
                if not chrom :
                    chrom = t[0][3]
                    idx = 10010
                elif chrom != t[0][3] :
                    chrom = t[0][3]
                    idx = 10010
                else :
                    idx += 10
                infor = "ID=AT" + pre + "-" + chrom + "G" + str(idx)
                if flag == "prot" :
                    infor = infor + ";Note=protein_coding_gene"
                else :
                    infor = infor + ";Note=transposable_element_gene"  
                mID = 0
                id = t[8].split(";")[0]
                id = id.split("=")[1] 
                newId =  "AT" + pre + "-" + chrom + "G" + str(idx)
                fo3.write("{}\t{}\n".format(id,newId))
                
            elif t[2] == "mRNA" :
                mID += 1
                infor = "ID=AT" + pre + "-" + chrom + "G" + str(idx) + "." + str(mID) + ";Parent=AT" + pre + "-" + chrom + "G" + str(idx)
                cdsID = 0
            elif t[2] =="exon" :
                nn = t[8].split(";")[0]
                nn = nn.split("=")[1]
                tt = nn.split(".")
                mm = tt[-1]
                cdsID = re.sub(r"exon","", mm)
                infor = "ID=AT" + pre + "-" + chrom + "G" + str(idx) + "." + str(mID) + "." + mm 
                infor = infor + ";Parent=AT" + pre + "-" + chrom + "G" + str(idx) + "." + str(mID)
            elif t[2] =="CDS" :
                nn = t[8].split(";")[0]
                nn = nn.split("=")[1]
                tt = nn.split(".")
                #cdsId = tt[-1]
                #infor = "ID=AT" + pre + "-" + chrom + "G" + str(idx) + "." + str(mID) +  ".cds" + str(cdsId) 
                infor = "ID=AT" + pre + "-" + chrom + "G" + str(idx) + "." + str(mID) +  ".cds" + str(cdsID)
                infor = infor + ";Parent=AT" + pre + "-" + chrom + "G" + str(idx) + "." + str(mID)                     
            outline ="\t".join(t[0:8]) + "\t" + infor + "\n"
            fo1.write(outline)  
            if flag == "prot" : fo2.write(outline)                        
    fi.close()
    fo1.close()
    fo2.close()
    fo3.close()
    
    
    ProtFas = os.path.join(outdir, acc + "." + ver + ".protein.fasta")
    CDSFas = os.path.join(outdir, acc + "." + ver + ".CDS.fasta")
    cDNAFas = os.path.join(outdir, acc + "." + ver + ".cDNA.fasta")
    GeneFas = os.path.join(outdir, acc + "." + ver + ".gene.fasta")
    
    EVM_Utils = "/projects/dep_coupland/grp_schneeberger/bin/EVM_r2012-06-25/EvmUtils/"
    cmd = EVM_Utils + "/gff3_file_to_proteins.pl " + outFile2  + " " + genome + " prot |awk '{if (/>/) print $1 ;else print }' > " + ProtFas
    print cmd
    os.system(cmd)
    cmd = EVM_Utils + "/gff3_file_to_proteins.pl " + outFile2  + " " + genome + " CDS  |awk '{if (/>/) print $1 ;else print }' > " + CDSFas
    os.system(cmd)
    cmd = EVM_Utils + "/gff3_file_to_proteins.pl " + outFile2  + " " + genome + " cDNA |awk '{if (/>/) print $1 ;else print }' > " + cDNAFas
    os.system(cmd)
    cmd = EVM_Utils + "/gff3_file_to_proteins.pl " + outFile2  + " " + genome + " gene |awk '{if (/>/) print $1 ;else print }' > " + GeneFas
    os.system(cmd)
    #os.system("rm " + outdir + "/genes.tmp.gff")
        
if __name__ == "__main__":
   main(sys.argv[1:])  
              
            
            
                                    
            