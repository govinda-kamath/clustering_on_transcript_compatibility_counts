
# coding: utf-8

# In[1]:

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
    opts, args = getopt.getopt(sys.argv[1:],"i:n:b:",["idir=","njobs=","base-dir="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python get_UMI_for_sampled_reads.py -i input_read_dir -b base-umi-dir [-n number-of-processes-to-use]')
    sys.exit(1)

read_dir=''
num_proc=1
UMI_dir=''
for opt,arg in opts:
    if opt in ("-i", "--idir"):
        read_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-b","--base-dir"):
        UMI_dir=arg

if (not read_dir) or (not UMI_dir) :
    print ('usage is : \n python get_UMI_for_sampled_reads.py -i input_read_dir [-n number-of-processes-to-use]')
    sys.exit(1)


# In[21]:

def get_umis(fltuple):
    read_dir=fltuple[0]
    flname=fltuple[1]
    cellname=flname.split('.')[0]
    #print cellname
    UMI_dir=fltuple[2]
    UMI_file=UMI_dir+cellname+'.umi'
    rnum_file=read_dir+flname
    out_file=read_dir+cellname+'.umi'
    
    cmd=""" awk '
    NR == FNR {
        for (i=1; i<=NF; i++) {
        linenums[$i]
            }
        }
    NR != FNR {
        if (FNR in linenums) {
        print
            }
        }
        ' """ + rnum_file +" "+ UMI_file+ " > "+ out_file
    #print cmd
    os.system(cmd)


# In[ ]:



fls=[x for x in os.listdir(read_dir) if x.endswith('.rnum')]
fltuple=itertools.product([read_dir],fls,[UMI_dir])
pool=mp.Pool(processes=num_proc)
pool.map(get_umis,fltuple)
#     for tup in fltuple:
#         get_umis(tup)