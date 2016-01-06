import os
import numpy as np
import sys
import itertools
import multiprocessing as mp
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:],"r:n:",["ref-transcriptome=","njobs="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python process_xprs.py  [-n number-of-processes-to-use]')
    sys.exit(1)

num_proc=1
for opt,arg in opts:
    if opt in ("-n","--njobs"):
        num_proc=int(arg)

def get_expression(fltuple):
    #print fltuple
    bowtie=fltuple[0]
    drname=fltuple[1]
    druse=bowtie+drname+'/'
    #print druse
    if 'SRR' in druse:
        flnames=os.listdir(druse)
        if 'results.xprs' in flnames:
            cmd = "cat "+druse+"results.xprs | awk '{print $2, $15}' | tail -n+2 | sort -k 1 > "+druse+"results.t3i" 
            os.system(cmd)

for samp in ['100','10','5','1','_point5','_point1']:
    bowtie='./Zeisel_Bowtie_subsample'+samp+'/'
    files = sorted(os.listdir(bowtie))
    #print files
    fltuple=itertools.product([bowtie],files)
    pool = mp.Pool(processes = num_proc)
    pool.map(get_expression,fltuple)

    
