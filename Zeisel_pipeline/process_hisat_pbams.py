
# coding: utf-8

# In[1]:

import os 
import itertools
import multiprocessing as mp
import numpy as np
import getopt
import sys


# In[3]:
if len(sys.argv) != 2:
    print("Usage is \n python process_hisat_pbams.py num_processes")

num_proc=int(sys.argv[1])
transcript_fl='./transcript_length.txt'


# In[4]:

trans_hash={}
index=0
with open(transcript_fl,'r') as f:
    for line in f:
        trans=line.split()[0]
        trans_hash[trans]=index
        index+=1


# In[21]:


pbam_dir='./hisat_pbams/'
out_dir='./hisat_t3i/'
os.system('mkdir -p '+out_dir)


# In[16]:

def process_pbam(fltuple):
    flname=fltuple[0]
    pbam_dir=fltuple[1]
    out_dir=fltuple[2]
    cellname=flname.split('.')[0]
    pbam_fl=pbam_dir+cellname+'.pbam'
    out_fl=out_dir+cellname+'.t3i'

    read_to_EQ={}
    with open(pbam_fl,'r') as f:
        for line in f:
            line1=line.split()
            read_id=int(line1[0])
            trans_id=line1[1]
            read_to_EQ.setdefault(read_id,set())
            read_to_EQ[read_id].add(trans_hash[trans_id])

    EQ_to_num_reads={}
    for read,EQ in read_to_EQ.items():
        tup_EQ=",".join(map(str,sorted(EQ)))
        EQ_to_num_reads.setdefault(tup_EQ,0)
        EQ_to_num_reads[tup_EQ]+=1

    with open(out_fl, 'w') as f:
        for EQ, vl in EQ_to_num_reads.items():
            f.write(EQ +"\t" +str(vl)+"\n")


# In[20]:

fl_list=os.listdir(pbam_dir)
fltuple=itertools.product(fl_list,[pbam_dir],[out_dir])

pool=mp.Pool(processes=num_proc)
pool.map(process_pbam,fltuple)


# In[ ]:



