# Run Bowtie on cells
import multiprocessing as mp
import numpy as np
import os
import sys
import getopt
import itertools


index=''
fldir=''
num_proc=1
out_base=''

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:o:n:r:",["idir=","odir=","njobs=","ref-index-path="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python run_bowtie.py -i input-read-dir -o output-kallisto-dir '+
           ' -r path-to-mouse-reference-index  [-n number-of-processes-to-use]')
    sys.exit(1)
    
    
for opt,arg in opts:
    if opt in ("-i", "--idir"):
        fldir=arg
    elif opt in ("-o","--odir"):
        out_base=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-r","--ref-index-path"):
        index=arg

if (not fldir) or (not index) or (not out_base):
    print ('usage is : \n python run_bowtie.py -i input-read-dir -o output-kallisto-dir '+
           ' -r path-to-mouse-reference-index  [-n number-of-processes-to-use]')
    sys.exit(1)

def run_bowtie1(fltuple):
    flname=fltuple[0]
    index = fltuple[1]
    fldir = fltuple[2]
    flpath = fldir+flname
    out_base = fltuple[3]
    out = out_base + flname.split('.')[0]
    #print (out)
    os.system('mkdir -p ' + out)
    BTcmd = 'gzip -dc '+flpath+' | bowtie -S -k 240 -m 1000 --offrate 1 '+index+' - | samtools view -Sb - > '+out+'/hits.bam'
    #print BTcmd
    os.system(BTcmd)

files = [x for x in os.listdir(fldir) if x.endswith('.fastq.gz')]
os.system('mkdir -p '+out_base)
#print files
fltuple=itertools.product(files,[index],[fldir],[out_base])
pool = mp.Pool(processes = num_proc)
pool.map(run_bowtie1,fltuple)

