# Run TopHat on each of selected cells
import multiprocessing as mp
import numpy as np
import itertools
import os
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:],"r:n:",["ref-transcriptome=","njobs="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python run_express.py -r ref-transcriptome  [-n number-of-processes-to-use]')
    sys.exit(1)
    
transcripts=''
num_proc=1
for opt,arg in opts:
    if opt in ("-r", "--ref-transcriptome"):
        transcripts=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
if (not transcripts):
    print ('usage is : \n python run_express.py -r ref-transcriptome  [-n number-of-processes-to-use]')
    sys.exit(1)

def run_express(fltuple):
    transcripts=fltuple[0]
    bam_dir=fltuple[1]
    flname=fltuple[2]
    #transcripts = '/data/SS_RNA_seq/Zeisel/reference_transcriptome/Mus_musculus.GRCm38.rel79.cdna.all.fa'
    hits = bam_dir+flname+'/'
    Ecmd = 'express --no-bias-correct -o '+hits+' ' + transcripts + ' ' + hits +'hits.bam'
    os.system(Ecmd)

def run_express100(fltuple)
    transcripts=fltuple[0]
    bam_dir=fltuple[1]
    flname=fltuple[2]
    #transcripts = '/data/SS_RNA_seq/Zeisel/reference_transcriptome/Mus_musculus.GRCm38.rel79.cdna.all.fa'
    hits = bam_dir+flname+'/'
    Ecmd = 'express -o '+hits+' ' + transcripts + ' ' + hits +'hits.bam'
    os.system(Ecmd)

for suffix in ['100','10','5','1','_point5','_point1']:
    bam_dir='./Zeisel_Bowtie_subsample'+suffix+'/'
    flnames=os.listdir(bam_dir)
    fltuple=itertools.product([transcripts],[bam_dir],flnames)
    if suffix=='100':
        pool = mp.Pool(processes = num_proc)
        pool.map(run_express100,fltuple)
    else:
        pool = mp.Pool(processes = num_proc)
        pool.map(run_express,fltuple)

