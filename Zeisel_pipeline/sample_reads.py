import numpy as np
import os
import commands
import itertools
import multiprocessing as mp
import scipy.sparse
import math
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:o:r:k:n:",["idir=","odir=","sampling-rate=","sampling-key=","njobs="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python sample_reads.py -i input_read_dir -o output-read-dir -r sampling-rate -k sampling-key [-n number-of-processes]')
    sys.exit(1)

key=0
read_dir=''
out_dir=''
sampling_rate=0
num_proc=1

print opts
for opt,arg in opts:
    if opt in ("-i", "--idir"):
        read_dir=arg
    elif opt in ("-o", "--odir"):
        out_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-r","--sampling-rate"):
        sampling_rate=float(arg)
    elif opt in ("-k","--sampling-key"):
        key=int(arg)
        
if (not read_dir) or (not out_dir) or (not sampling_rate) or (not key):
    print ('usage is : \n python sample_reads.py -i input_read_dir -o output-read-dir -r sampling-rate -k sampling-key [-n number-of-processes]')
    sys.exit(1)
    


def subsample_reads(fltuple):
    flname=fltuple[0]
    subsampling_rate=fltuple[1]
    read_dir=fltuple[2]
    out_dir=fltuple[3]
    key=fltuple[4]
    
    cmd='zcat '+read_dir+flname+' | wc -l'
    out=commands.getstatusoutput(cmd)
    
    num_reads=int(out[1])/4
    num_reads_to_sample=int(math.ceil(num_reads*subsampling_rate))
    print num_reads_to_sample
    
    cmd2='seqtk sample -s'+str(key)+' '+read_dir+flname+' '+ str(num_reads_to_sample) + ' | gzip > '+out_dir+flname
    print cmd2
    os.system(cmd2)


# In[14]:

#print repr(sys.argv[1])
os.system('mkdir -p '+out_dir)
flnames=sorted([x for x in os.listdir(read_dir) if x.endswith('.fastq.gz')])
fltuple=itertools.product(flnames,[sampling_rate],[read_dir],[out_dir],[key])
pool=mp.Pool(processes=num_proc)
pool.map(subsample_reads,fltuple)

# In[15]: