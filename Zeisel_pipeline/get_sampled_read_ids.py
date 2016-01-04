# coding: utf-8

# In[1]:


import os
import itertools
import multiprocessing as mp
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:",["idir=","njobs="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python get_sampled_read_ids.py -i input_read_dir [-n number-of-processes-to-use]')
    sys.exit(1)

read_dir=''
num_proc=1
for opt,arg in opts:
    if opt in ("-i", "--idir"):
        read_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)

if not read_dir:
    print ('usage is : \n python get_sampled_read_ids.py -i input_read_dir [-n number-of-processes-to-use]')
    sys.exit(1)

# In[4]:

def get_read_ids(fltuple):
    read_dir=fltuple[1]
    flname=fltuple[0]
    cellname=flname.split('.')[0]
    cmd="zcat "+read_dir+flname + " | awk -F'[ .]' 'NR%4==1 {print $2}' > "+read_dir+cellname+".rnum"
    #print cellname
    os.system(cmd)
    


# In[5]:

    
flnames=sorted([x for x in os.listdir(read_dir) if x.endswith('.fastq.gz')])
fltuple=itertools.product(flnames,[read_dir])
pool=mp.Pool(processes=64)
pool.map(get_read_ids,fltuple)


# In[ ]:



