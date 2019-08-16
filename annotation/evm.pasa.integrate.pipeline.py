#!/usr/bin/env python
# encoding: utf-8
'''
Created on Jul 26, 2017

@author: jiao
'''
import sys
import os
import getopt
import re
import glob
import subprocess
import time
from util.util import *


def main(argv):
    cfgFile = ""
    try:
        opts, args = getopt.getopt(argv,"f:",["cfg=",])  ##outfile: AGP file, Bed file, chr-fasta file
    except getopt.GetoptError:
        print 'evm.pasa.integrate.pipeline.py -f <cfg>  '
        sys.exit(2)
    if len(opts) == 0 :
        print 'evm.pasa.integrate.pipeline.py -f <cfg>  '
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print 'evm.pasa.integrate.pipeline.py -f <cfg> '
            sys.exit()
        elif opt in ("-f", "--cfg"):
            cfgFile = arg

    runCfg = parserConfig(cfgFile)
    RUNLOG = os.path.join( runCfg["workDir"], "run.log" )
    LOG = open(RUNLOG,"w", 1)

    ''' step 1: run exonerate to align protein sequences '''
    protAlnOutdir = os.path.join( runCfg["protOutdir"], runCfg["protAligner"] )
    splitProt(runCfg["protDir"],runCfg["protChunkNum"],protAlnOutdir)

    protFas = glob.glob(protAlnOutdir + "/*.fa")
    protFas.sort()
    bjobs = []
    exoOutFiles = []
    LOG.write("------ step1: run exonerate  ------\n")
    for i in range(len(protFas)) :
        #break
        name = os.path.basename(protFas[i])
        name = name.split(".")[2]
        cmdExonerate = "exonerate --ryo \'AveragePercentIdentity: %pi\\n\' --showvulgar no --showalignment no --showquerygff no --model protein2genome --showtargetgff yes "
        cmdExonerate = cmdExonerate + " --percent " + runCfg["exoP"] +  " --minintron " + runCfg["exoMinIntron"] + " --maxintron " + runCfg["exoMaxIntron"]
        exoOut = os.path.join(protAlnOutdir, "ex.chunk." + str(name) + ".out")
        bsubout = os.path.join(protAlnOutdir, "bsub_out." + str(name) + ".log")
        bsuberr = os.path.join(protAlnOutdir, "bsub_err." + str(name) + ".log")
        bsub = "bsub -q normal -n 1 -R'rusage[mem=2000]' "
        bsub = bsub + " -o " + bsubout + " -e " + bsuberr
        cmdExonerate = cmdExonerate + " " + protFas[i] + " " + runCfg["genome"] + " > " + exoOut
        runJob = bsub + " \" " + cmdExonerate + " \" "
        LOG.write(runJob + "\n")
        x = subprocess.check_output(runJob,shell=True)
        #mat = re.search(r'<(\d+)>',x)
        #jobID = mat.group(1)
        #bjobs.append(jobID)
        bjobs.append(bsubout)
        exoOutFiles.append(exoOut)
    LOG.write("\n======\n\n")

    ''' step 2; run hist2, stringTie to map RNA-seq reads and assembly transcripts '''
    LOG.write("------ step2: run hisat2, stringTie to map RNA-seq reads and assemble transcripts ------\n")
    rnaseqFastqs1 = glob.glob(runCfg["rnaseqDataDir"] + "/*_1.fastq*")
    rnaseqFastqs2 = glob.glob(runCfg["rnaseqDataDir"] + "/*_2.fastq*")
    rnaseqFastqs0 = glob.glob(runCfg["rnaseqDataDir"] + "/*_trim.fastq*")
    #rnaseqFastqs0 = filter(lambda x: '_' not in x , rnaseqFastqs0)
    if len(rnaseqFastqs1) > 0 :
        rnaseqFastqs1.sort()
        rnaseqFastqs2.sort()
    if len(rnaseqFastqs0) > 0 :
        rnaseqFastqs0.sort()
    if len(rnaseqFastqs0)==0 and len(rnaseqFastqs1)==0 :
        LOG.write("no rnaseq data \n")
    merlist =os.path.join(runCfg["rnaseqOutdir"],  "mergelist.txt")
    if (os.path.isfile(merlist)) : os.system("rm " + merlist)
    fo = open(merlist,"w")
    rnaseqFastqs = rnaseqFastqs1 + rnaseqFastqs0
    for i in range(len(rnaseqFastqs)) :
        #break
        n = i + 1
        read1 = rnaseqFastqs[i]
        if (len(rnaseqFastqs2) > i ) :
            read2 = rnaseqFastqs2[i]
        else :
            read2 = "None"
        print read1, read2
        name = os.path.basename(read1)
        name = name.split("_")[0]
        jobShell = os.path.join(runCfg["rnaseqOutdir"],"run." + name + ".sh")
        bsubout = os.path.join(runCfg["rnaseqOutdir"], "bsub_out." + name + ".log")
        bsuberr = os.path.join(runCfg["rnaseqOutdir"], "bsub_err." + name + ".log")
        bsub = "bsub -q multicore20 -R'rusage[mem=32000]span[hosts=1]' -n " + runCfg["rnaThread"]
        if int(runCfg["rnaThread"]) > 20 :
            bsub = "bsub -q multicore40 -R'rusage[mem=32000]span[hosts=1]' -n " + runCfg["rnaThread"]
        bsub = bsub + " -o " + bsubout + " -e " + bsuberr
        runHisatStringtie (name,read1,read2,runCfg["hisat2index"],runCfg["rnaseqOutdir"],runCfg["rnaThread"],jobShell)
        runJob = bsub + " bash " + jobShell
        fo.write(os.path.join(runCfg["rnaseqOutdir"], name + ".gtf") + "\n")
        LOG.write(runJob + "\n")
        x = subprocess.check_output(runJob,shell=True)
        #mat = re.search(r'<(\d+)>',x)
        #jobID = mat.group(1)
        #bjobs.append(jobID)
        bjobs.append(bsubout)
    LOG.write("\n======\n\n")
    fo.close()


    ''' step 3: run augustus,glimmerhmm, snap ab initio '''
    LOG.write("------step3: run augustus,glimmerhmm/genemark.hmm ab initio ------\n")
    splitGenome(runCfg["genome"],runCfg["genomeCSize"],runCfg["genSpOutdir"])
    genomeFas = glob.glob(runCfg["genSpOutdir"]+"/genome.chunk.*fa")
    genomeFas.sort()
    augOut1 = []
    gliOut = []
    snaOut = []

    for i in range(len(genomeFas)) :
        #break
        name = os.path.basename(genomeFas[i])
        name = name.split(".")[2]
        bsubout = os.path.join(runCfg["augOutdir"], "bsub_out." + str(name) + ".log")
        bsuberr = os.path.join(runCfg["augOutdir"], "bsub_err." + str(name) + ".log")
        bsub = "bsub -q normal -n 1 -R'rusage[mem=2000]' -o " + bsubout + " -e " + bsuberr
        outFile = os.path.join(runCfg["augOutdir"], "augustus.no-hint.chunk." + str(name) + ".gff")
        augCmd = "augustus --UTR=on --print_utr=on --exonnames=on --codingseq=on  --genemodel=complete --alternatives-from-evidence=true --gff3=on"
        augCmd = augCmd + " --species=" + runCfg["augMSpecies"] + " --extrinsicCfgFile=" + runCfg["augAbConfig"] + " --outfile=" + outFile + " " + genomeFas[i]
        runJob = bsub + " " +  augCmd
        print "run job: ", runJob
        LOG.write(runJob + "\n")
        x = subprocess.check_output(runJob,shell=True)
        #mat = re.search(r'<(\d+)>',x)
        #jobID = mat.group(1)
        augOut1.append(outFile)
        #bjobs.append(jobID)
        bjobs.append(bsubout)

        #run glimmerhmm
        bsubout = os.path.join(runCfg["gliOutdir"], "bsub_out." + str(name) + ".log")
        bsuberr = os.path.join(runCfg["gliOutdir"], "bsub_err." + str(name) + ".log")
        bsub = "bsub -q normal -n 1 -R'rusage[mem=2000]' -o " + bsubout + " -e " + bsuberr
        outFile = os.path.join(runCfg["gliOutdir"], "glimmer.chunk." + str(name) + ".gff")
        gliCmd = "glimmerhmm " + genomeFas[i] + " -d " + runCfg["gliTrainDir"] + " -o " + outFile +  " -g -f"
        runJob = bsub + " " +  gliCmd
        print "run job: ", runJob
        LOG.write(runJob + "\n")
        x = subprocess.check_output(runJob,shell=True)
        #mat = re.search(r'<(\d+)>',x)
        #jobID = mat.group(1)
        gliOut.append(outFile)
        #bjobs.append(jobID)
        bjobs.append(bsubout)

        #run SNAP
        bsubout = os.path.join(runCfg["SNAPOutdir"], "bsub_out." + str(name) + ".log")
        bsuberr = os.path.join(runCfg["SNAPOutdir"], "bsub_err." + str(name) + ".log")
        bsub = "bsub -q normal -n 1 -R'rusage[mem=2000]' -o " + bsubout + " -e " + bsuberr
        outFile = os.path.join(runCfg["SNAPOutdir"], "SNAP.chunk." + str(name) + ".gff")
        snapCmd = "snap " + runCfg["SNAPhmm"] + " " + genomeFas[i] + " -gff > " + outFile
        runJob = bsub + " "  + "\"" + snapCmd + "\""
        print "run job: ", runJob
        LOG.write(runJob + "\n")
        x = subprocess.check_output(runJob,shell=True)
        #mat = re.search(r'<(\d+)>',x)
        #jobID = mat.group(1)
        snaOut.append(outFile)
        #bjobs.append(jobID)
        bjobs.append(bsubout)

        #run genemarkHMM

    LOG.write("\n======\n\n")
    jobWait(bjobs)
    bjobs = []

    ''' step 4: merge result of each kind of evidence '''
    LOG.write("------ step 4: merge result of each kind of evidence ------\n")
    ##1) generate hint from exonerate outputs for augustus
    exoMergeOut = os.path.join(runCfg["protOutdir"],"exonerate.out")
    if os.path.isfile(exoMergeOut) : os.system("rm " + exoMergeOut)

    for exo in exoOutFiles:
        break
        os.system("cat " + exo + " >> " + exoMergeOut)
    exohint = os.path.join(runCfg["protOutdir"],"exonerate.hints")
    cmd = "exonerate2hints.pl --in=" + exoMergeOut + " --source=P --out=" + exohint + " --minintronlen=" + runCfg["exoHintMinIntron"] + " --maxintronlen=" + runCfg["exoHintMaxIntron"]
    print "run ... ", cmd, "\n"
    LOG.write(cmd + "\n")
    os.system(cmd )

    ##generate protein-evidence gff for evm
    exo4evmPerl = runCfg["EVM_misc"] + "/exonerate_gff_to_alignment_gff3.pl"
    exo4evmGff = os.path.join(runCfg["protOutdir"],"exonerate.4evm.gff")
    cmd = exo4evmPerl + " " + exoMergeOut + " > " + exo4evmGff
    print "run ... ", cmd, "\n"
    LOG.write(cmd + "\n")
    os.system(cmd)

    ##2) generat hint from stringtie gtf files
    rnalist =os.path.join(runCfg["rnaseqOutdir"], "mergelist.txt")
    rnagtf = os.path.join(runCfg["rnaseqOutdir"], "stringtie.merged.gtf")
    rnahint = os.path.join(runCfg["rnaseqOutdir"], "rnaseq.stringtie.hints")
    if len(rnaseqFastqs) >0 :
        os.system("stringtie --merge -p 1 -o " + rnagtf + " " + rnalist)
        striGTF2augHint(rnagtf,rnahint)
    ##generate RNAseq-evidenc gff for evm
        rna4evmPerl = runCfg["TransDeUtil"] + "/cufflinks_gtf_to_alignment_gff3.pl"
        rna4evmGff = os.path.join(runCfg["rnaseqOutdir"], "rnaseq.4evm.gff")
        cmd =rna4evmPerl + " " + rnagtf + " > " + rna4evmGff
        LOG.write(cmd + "\n")
        os.system(cmd)

    ##3) transfer the gff output of ab ininitio prediction into gff for evm
    augMergeOut = os.path.join(runCfg["abiOutdir"],"augustus.merge.out")
    gliMergeOut = os.path.join(runCfg["abiOutdir"],"glimmerhmm.merge.out")
    snaMergeOut = os.path.join(runCfg["abiOutdir"],"snap.merge.out")
    if os.path.isfile(augMergeOut) : os.system("rm "+ augMergeOut)
    if os.path.isfile(gliMergeOut) : os.system("rm "+ gliMergeOut)
    if os.path.isfile(snaMergeOut) : os.system("rm "+ snaMergeOut)
    for i in range(len(genomeFas)) :
        break
        os.system("cat " + augOut1[i] + " >> " + augMergeOut)
        os.system("cat " + gliOut[i] + " >> " + gliMergeOut)
        os.system("cat " + snaOut[i] + " >> " + snaMergeOut)

    aug4evmPerl = runCfg["EVM_misc"] + "/augustus_to_GFF3.pl"
    aug4evmGff = os.path.join(runCfg["abiOutdir"], "augustus.4evm.gff")
    cmd = aug4evmPerl + " " + augMergeOut + " > " + aug4evmGff
    LOG.write(cmd + "\n")
    os.system(cmd)

    gli4evmPerl = runCfg["EVM_misc"] + "/glimmerHMM_to_GFF3.pl"
    gli4evmGff = os.path.join(runCfg["abiOutdir"], "glimmerhmm.4evm.gff")
    cmd = gli4evmPerl + " " + gliMergeOut + " > " + gli4evmGff
    LOG.write(cmd + "\n")
    os.system(cmd)

    sna4evmPerl = runCfg["EVM_misc"] + "/SNAP_to_GFF3.jiao.pl"
    sna4evmGff = os.path.join(runCfg["abiOutdir"], "SNAP.4evm.gff")
    cmd = sna4evmPerl + " " + snaMergeOut + " > " + sna4evmGff
    LOG.write(cmd + "\n")
    os.system(cmd)

    abi4evmGff = os.path.join(runCfg["abiOutdir"], "abinitio.4evm.gff")
    LOG.write(cmd + "\n")
    os.system("cat " + aug4evmGff + " " + gli4evmGff + " " + sna4evmGff + " > " + abi4evmGff )

    LOG.write("\n======\n\n")

    ''' step 5: run augustus with evidence '''
    LOG.write("------ step 5: run augustus with evidence ------\n")
    bjobID = []
    bjobs = []
    for i in range(len(genomeFas)) :
        #break
        n = i + 1
        bsubout = os.path.join(runCfg["augOutdir"], "bsub_out." + str(n) + ".hint.log")
        bsuberr = os.path.join(runCfg["augOutdir"], "bsub_err." + str(n) + ".hint.log")
        bsub = "bsub -q normal -n 1 -R'rusage[mem=2000]' "
        bsub = bsub + " -o " + bsubout + " -e " + bsuberr
        outFile = os.path.join(runCfg["augOutdir"], "augustus.hint.chunk." + str(n) + ".gff")
        augCmd = "augustus --UTR=on --print_utr=on --exonnames=on --codingseq=on  --genemodel=complete --alternatives-from-evidence=true --gff3=on"
        augCmd = augCmd + " --species=" + runCfg["augMSpecies"] + " --extrinsicCfgFile=" + runCfg["augHintConfig"] + " --outfile=" + outFile + " " + genomeFas[i]
        runJob = bsub + " " +  augCmd
        LOG.write(runJob + "\n")
        x = subprocess.check_output(runJob,shell=True)
        #mat = re.search(r'<(\d+)>',x)
        #jobID = mat.group(1)
        #bjobID.append(jobID)
        bjobs.append(bsubout)
    LOG.write("\n======\n\n")
    jobWait(bjobs)
    bjobID = []
    bjobs = []

    ''' step 5: run EVM + PASA '''
    LOG.write("------ step 5: run EVM + PASA ------\n")
    EVM_Utils = runCfg["EVM_Utils"]

    ##Partitioning the Inputs
    evmSplitdir = runCfg["evmOutdir"] + "/split"
    os.system("mkdir -p " + evmSplitdir)
    #else :
    #    os.system("rm -r " + evmSplitdir + "/*")

    os.chdir(evmSplitdir)
    partlistOut = os.path.join(evmSplitdir, "evm.partition.list")
    cmd = EVM_Utils + "/partition_EVM_inputs.pl --segmentSize 100000 --overlapSize 10000  --genome " + runCfg["genome"]
    cmdPart =  " --gene_predictions " + abi4evmGff
    cmdPart = cmdPart + " --protein_alignments " + exo4evmGff
    if len(rnaseqFastqs) > 0 :
        cmdPart = cmdPart + " --transcript_alignment " + rna4evmGff
    cmd = cmd + " " + cmdPart + " --partition_listing " + partlistOut
    LOG.write(cmd + "\n")
    os.system(cmd)

    ##Generating the EVM Command Set
    cmdlist = os.path.join(evmSplitdir, "evm.commands.list")
    evmOut = "evm.out"
    cmd = EVM_Utils + "/write_EVM_commands.pl "
    cmd = cmd + " --genome " + runCfg["genome"] + " --weights " + runCfg["evmWeight"]
    cmd = cmd + " " + cmdPart + " --partitions " + partlistOut + " --output_file_name " + evmOut + " > " + cmdlist
    LOG.write(cmd + "\n")
    os.system(cmd)

    bjobs = []
    fi = open(cmdlist, "r")
    chunkOut = os.path.join(evmSplitdir, "evm.cmd.chunk.0.sh")
    fo = open( chunkOut, "w")

    cmdNum = int(runCfg["evmCmdNum"])
    k = 0
    n = 0
    while True :
        line = fi.readline()
        if not line : break
        if k == cmdNum :
            fo.close()

            bsubout = os.path.join(evmSplitdir , "bsub_out." + str(n) + ".log")
            bsuberr = os.path.join(evmSplitdir , "bsub_err." + str(n) + ".log")
            if os.path.isfile(bsubout) :
                os.system("rm -r " + bsubout)
                os.system("rm -r " + bsuberr)
            bsub = "bsub -q normal -n 1 -R'rusage[mem=8000]' "
            bsub = bsub + " -o " + bsubout + " -e " + bsuberr
            runJob = bsub + " bash " + chunkOut
            LOG.write(runJob + "\n")
            x = subprocess.check_output(runJob,shell=True)
            #mat = re.search(r'<(\d+)>',x)
            #jobID = mat.group(1)
            #bjobID.append(jobID)
            bjobs.append(bsubout)
            n = n + 1
            chunkOut = os.path.join(evmSplitdir, "evm.cmd.chunk." + str(n) + ".sh")
            fo = open( chunkOut, "w")
            fo.write(line)
            k = 1

        else :
            fo.write(line)
            k += 1
    fi.close()
    if k > 0 :
        fo.close()
        bsubout = os.path.join(evmSplitdir, "bsub_out." + str(n) + ".log")
        bsuberr = os.path.join(evmSplitdir, "bsub_err." + str(n) + ".log")
        bsub = "bsub -q normal -n 1 -R'rusage[mem=8000]' "
        bsub = bsub + " -o " + bsubout + " -e " + bsuberr
        runJob = bsub + " bash " + chunkOut
        LOG.write(runJob + "\n")
        x = subprocess.check_output(runJob,shell=True)
        #mat = re.search(r'<(\d+)>',x)
        #jobID = mat.group(1)
        #bjobID.append(jobID)
        bjobs.append(bsubout)

    jobWait(bjobs)
    bjobs = []

    ##Combining the Partitions
    cmd = EVM_Utils + "/recombine_EVM_partial_outputs.pl  --partitions " + partlistOut + " --output_file_name evm.out"
    LOG.write(cmd + "\n")
    os.system(cmd)
    ##Convert to GFF3 Format
    cmd = EVM_Utils + "/convert_EVM_outputs_to_GFF3.pl --partitions " + partlistOut + " --output evm.out  --genome " + runCfg["genome"]
    LOG.write(cmd + "\n")
    os.system(cmd)
    ##merge the out files and gff files
    evmOutFiles1 = glob.glob(evmSplitdir + "/chr*/evm.out")
    evmOutFiles1.sort()
    evmOutFiles2 = glob.glob(evmSplitdir + "/*/evm.out")
    evmOutFiles2 = filter(lambda x: 'chr' not in x , evmOutFiles2)
    evmOutFiles2.sort()
    evmOutFiles = evmOutFiles1 + evmOutFiles2

    evmGffFiles1 = glob.glob(evmSplitdir + "/chr*/evm.out.gff3")
    evmGffFiles1.sort()
    evmGffFiles2 = glob.glob(evmSplitdir + "/*/evm.out.gff3")
    evmGffFiles2 = filter(lambda x: 'chr' not in x , evmGffFiles2)
    evmGffFiles2.sort()
    evmGffFiles = evmGffFiles1 + evmGffFiles2

    evmOutAll = os.path.join(runCfg["evmOutdir"], "evm.all.out")
    evmGffAll = os.path.join(runCfg["evmOutdir"], "evm.all.gff3")
    if os.path.isfile(evmOutAll) : os.system("rm " + evmOutAll)
    if os.path.isfile(evmGffAll) : os.system("rm " + evmGffAll)
    for i in range(len(evmOutFiles)) :
        os.system("cat " + evmOutFiles[i] + " >> " + evmOutAll)
        os.system("cat " + evmGffFiles[i] + " >> " + evmGffAll)

    evmProtFas = os.path.join(runCfg["evmOutdir"], "evm.annotation.protein.fasta")
    evmCDSFas = os.path.join(runCfg["evmOutdir"], "evm.annotation.CDS.fasta")
    evmcDNAFas = os.path.join(runCfg["evmOutdir"], "evm.annotation.cDNA.fasta")
    evmGeneFas = os.path.join(runCfg["evmOutdir"], "evm.annotation.gene.fasta")

    cmd = EVM_Utils + "/gff3_file_to_proteins.pl " + evmGffAll  + " " + runCfg["genome"] + " prot > " + evmProtFas
    LOG.write(cmd + "\n")
    os.system(cmd)
    cmd = EVM_Utils + "/gff3_file_to_proteins.pl " + evmGffAll  + " " + runCfg["genome"] + " CDS > " + evmCDSFas
    LOG.write(cmd + "\n")
    os.system(cmd)
    cmd = EVM_Utils + "/gff3_file_to_proteins.pl " + evmGffAll  + " " + runCfg["genome"] + " cDNA > " + evmcDNAFas
    LOG.write(cmd + "\n")
    os.system(cmd)
    cmd = EVM_Utils + "/gff3_file_to_proteins.pl " + evmGffAll  + " " + runCfg["genome"] + " gene > " + evmGeneFas
    LOG.write(cmd + "\n")
    os.system(cmd)

    ## use PASA to update the EVM consensus predictions, adding UTR annotations and models for alternatively spliced isoforms (leveraging D and E).
    #Two rounds

    LOG.write("\n======\n\n")
    LOG.write("Done\n")
    LOG.close()



if __name__ == "__main__":
   main(sys.argv[1:])


