#!/usr/bin/env python
# encoding: utf-8


import sys
import os
import getopt
import re
import glob
import multiprocessing


def main(argv):
    indir = ""
    outdir = ""
    ref = ""
    try:
        opts, args = getopt.getopt(argv,"i:o:r:",["indir=","outdir","ref="]) 
    except getopt.GetoptError:
        print 'run.scipio.py -i <indir> -o <outdir> -r <ref>'
        sys.exit(2)
    if len(opts) == 0 :
        print 'run.scipio.py -i <indir> -o <outdir> -r <ref>'
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print 'run.scipio.py -i <indir> -o <outdir> -r <ref>'
            sys.exit()
        elif opt in ("-i", "--indir"):
            indir = arg           
        elif opt in ("-o", "--outdir"):
            outdir = arg
        elif opt in ("-r", "--ref"):
            ref = arg
            
## gene model from augustus protein-based, scipio protein-based or stringTie-gff based
    if not os.path.isdir(outdir) :
        os.mkdir(outdir)    
    protFas = glob.glob(indir + "/*.fa")
    protFas.sort()
    pool = multiprocessing.Pool(processes = 60)
    outGff2 = []
    for file in protFas :
        name = os.path.basename(file)
        #print file
        out1 = outdir + "/" + name + ".scipio.out"
        out2 = outdir + "/" + name + ".scipio.gff"                
        pool.apply_async(runSci, (ref,file,out1,out2,))
        outGff2.append(out2)
    pool.close()
    pool.join()
        
def runSci(ref,prot,out1,out2):       
    sci = "/srv/netscratch/dep_coupland/grp_schneeberger/bin/scipio/scipio-1.4/scipio.1.4.1.pl"
    yam = "/srv/netscratch/dep_coupland/grp_schneeberger/bin/scipio/scipio-1.4/yaml2gff.1.4.pl"
    print  sci + " " + ref + " " + prot + " > " + out1
    os.system(sci + " " + ref + " " + prot + " > " + out1)
    os.system(yam + " " + out1 + " > " + out2)
    
    
    
    

if __name__ == "__main__":
   main(sys.argv[1:])     