
# coding: utf-8

# In[1]:

import os 
import itertools
import multiprocessing as mp
import numpy as np
import getopt
import sys

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:u:o:",["idir=","njobs=","UMI=","odir="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python process_pbam.py -i input-pbam-dir -o output-dir -u base-umi-dir [-n number-of-processes-to-use]')
    sys.exit(1)
    
UMI_dir = ''
pbam_dir=''
transcript_fl='./transcript_length.txt'
out_dir=''
num_proc=1

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        pbam_dir=arg
    elif opt in ("-o", "--odir"):
        out_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-u","--UMI"):
        UMI_dir=arg

if (not UMI_dir) or (not pbam_dir) or (not out_dir):
    print ('usage is : \n python process_pbam.py -i input-pbam-dir -o output-dir -u base-umi-dir [-n number-of-processes-to-use]')
    sys.exit(1)

# In[2]:

def make_trans_hash_global():
    global trans_hash


# In[3]:

def process_pbam(fltuple):
    flname=fltuple[0]
    UMI_dir=fltuple[1]
    pbam_dir=fltuple[2]
    out_dir=fltuple[3]
    
    cellname=flname.split('.')[0]
    pbam_fl=pbam_dir+cellname+'.pbam'
    UMI_fl=UMI_dir+cellname+'.umi'
    out_fl=out_dir+cellname+'.t3i'
    
    print cellname
    
    read_to_UMI=np.loadtxt(UMI_fl,usecols=(0,), dtype=str)
    
    read_to_EQ={}
    with open(pbam_fl,'r') as f:
        for line in f:
            line1=line.split()
            read_id=int(line1[0])
            trans_id=line1[1]
            read_to_EQ.setdefault(read_id,set())
            read_to_EQ[read_id].add(trans_hash[trans_id])
    
    trans_to_UMI={}
    for read_id,EQ in read_to_EQ.items():
        tup_EQ=",".join(map(str,sorted(EQ)))
        trans_to_UMI.setdefault(tup_EQ,set())
        trans_to_UMI[tup_EQ].add(read_to_UMI[read_id-1])
        
    EQ_to_num={}
    for EQ,UMI_set in trans_to_UMI.items():
        EQ_to_num[EQ]=len(UMI_set)
    
    with open(out_fl, 'w') as f:
        for EQ, vl in EQ_to_num.items():
            f.write(EQ +"\t" +str(vl)+"\n")



trans_hash={}
index=0
with open(transcript_fl,'r') as f:
    for line in f:
        trans=line.split()[0]
        trans_hash[trans]=index
        index+=1

make_trans_hash_global()


# In[6]:


os.system('mkdir -p '+out_dir)
flnames=sorted([x for x in os.listdir(pbam_dir) if x.endswith('.pbam')])
fltuple=itertools.product(flnames,[UMI_dir],[pbam_dir],[out_dir])
pool=mp.Pool(processes=num_proc)
pool.map(process_pbam,fltuple)
