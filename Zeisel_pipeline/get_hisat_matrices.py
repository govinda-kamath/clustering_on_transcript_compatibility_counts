# coding: utf-8

# In[3]:

import os 
import multiprocessing as mp
import numpy as np
import pickle
import scipy.sparse
import getopt
import sys

t3i_dir=''
out_file=''
out_ufile=''

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:t:d:",["idir=","TCC-file=","TCCD-file="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python get_hisat_matices.py -i input_t3i_dir  -t path-to-output-TCC-UMI-file -d path-to-output-TCC-UMI-dist-file')
    sys.exit(1)

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        t3i_dir=arg
    elif opt in ("-t","--TCC-file"):
        out_ufile=arg
    elif opt in ("-d","--TCCD-file"):
        out_file=arg
        
if (not t3i_dir) or (not out_ufile) or (not out_file):
    print ('usage is : \n python get_hisat_matices.py -i input_t3i_dir  -t path-to-output-TCC-UMI-file -d path-to-output-TCC-UMI-dist-file')
    sys.exit(1)

# In[4]:

def data_to_dist(X):
    s = np.sum(X,axis=1)
    X = X / s[:,None]
    return X, s



# In[5]:

def get_UMI_matrix(t3i_dir,out_file,out_ufile):
    
  

    flnames=sorted([x for x in os.listdir(t3i_dir) if x.endswith(".t3i")])

    eq_hash={}
    index=0
    for flname in flnames:
        with open(t3i_dir+flname) as f:
            for line in f:
                transid=line.split()[0]
                if transid not in eq_hash:
                    eq_hash[transid]=index
                    index+=1

    eq_mat=np.zeros((len(flnames),len(eq_hash)))
    for  index in xrange(len(flnames)):
        with open(t3i_dir+flnames[index]) as f:
            #print index
            for line in f:
                transid,num_obs=line.split()[0:2]
                index2=eq_hash[transid]
                eq_mat[index,index2] = int(num_obs)

    eq_mat_norm, num_mapped_reads =data_to_dist(eq_mat)
    S=scipy.sparse.csr_matrix(eq_mat_norm)
    Su=scipy.sparse.csr_matrix(eq_mat)


    with open(out_file, 'wb') as outfile:
        pickle.dump(S, outfile, pickle.HIGHEST_PROTOCOL)

    with open(out_ufile, 'wb') as outfile:
        pickle.dump(Su, outfile, pickle.HIGHEST_PROTOCOL)


# In[ ]:
get_UMI_matrix(t3i_dir,out_file,out_ufile)
