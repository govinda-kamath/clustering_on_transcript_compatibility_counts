
# coding: utf-8

# In[47]:

import os
import itertools
import multiprocessing as mp
import sys

if len(sys.argv) !=2:
    print('Usage is \n python process_hisat_sam.py num_procs')

# In[48]:

ip_dir='./hisat/'
op_dir='./hisat_pbams/'

os.system('mkdir -p '+op_dir)
num_proc=int(sys.argv[1])
samfiles=[flnm.split('.')[0] for flnm in os.listdir(ip_dir) if flnm.endswith('.sam')]


# In[49]:

def process_hisat_samfiles(fltuple):
    flname=fltuple[0]
    ippath=fltuple[1]
    oppath=fltuple[2]
    ip_file_path=ippath+flname+'.sam'
    op_file_path=oppath+flname+'.pbam'
    cmd=("grep -v ^@ "+ip_file_path+" | awk '{split($1,a,"+' "."); '
    + "print a[2], $3}' | grep -v '*' > "
    + op_file_path)
    os.system(cmd)


# In[50]:

fltuple=itertools.product(samfiles,[ip_dir], [op_dir])
pool=pool = mp.Pool(processes = num_proc)
pool.map(process_hisat_samfiles,fltuple)


# In[ ]:



