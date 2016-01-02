
# coding: utf-8

# In[1]:

import os 
import subprocess
import itertools
import multiprocessing as mp
import getopt
import sys


try:
    opts, args = getopt.getopt(sys.argv[1:],"i:o:n:",["idir=","odir=","njobs="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python process_SRA.py -i input_SRA_dir -o output_fastq_dir [-n number-of-processes-to-use]')
    sys.exit(1)


SRA_dir=''
out_dir=''
num_proc=1
for opt,arg in opts:
    #print (opt)
    #print (arg)
    if opt in ("-i", "--idir"):
        SRA_dir=arg
    elif opt in ("-o","--odir"):
        out_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)

if (not SRA_dir) or (not out_dir):
    print ('usage is : \n python process_SRA.py -i input_SRA_dir -o output_fastq_dir [-n number-of-processes-to-use]')
    sys.exit(1)

#print (SRA_dir)
#print (out_dir)
print("using "+str(num_proc)+" processes.")
#SRA_dir= '/data/SS_RNA_seq/Trapnell/aspera_DL/'
#out_dir='/data/SS_RNA_seq/Trapnell/read_data_test/'


# In[ ]:




# In[16]:

def process_SRA(fltuple):
    flname=fltuple[0]
    out_dir=fltuple[1]
    dir_path=fltuple[2]   
    
    cmd='fastq-dump --split-files --gzip '+ dir_path+flname+ ' -O '+ out_dir
    
    os.system(cmd)
    #print(cmd)


# In[15]:

flnames=sorted(os.listdir(SRA_dir))

fltuple=itertools.product(flnames, [out_dir],[SRA_dir])
os.system('mkdir -p '+out_dir)


# In[ ]:

pool=mp.Pool(processes=num_proc)
pool.map(process_SRA, fltuple)


# In[ ]:



