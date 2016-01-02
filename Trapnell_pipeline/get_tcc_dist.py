
# coding: utf-8

# In[1]:

import numpy as np
import os
import sys
import multiprocessing
import scipy.sparse
import pickle
import getopt


# In[2]:

def data_to_dist(X):
    s = np.sum(X,axis=1)
    X = X / s[:,None]
    return X, s


# In[3]:

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:m:t:d:",["idir=","num-eq-classes=","TCC-file=","TCCD-file="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python get_tcc_dist.py -i input_tcc_dir -m number-of-eq-classes -t path-to-output-TCC-file -d path-to-output-TCC-dist-file')
    sys.exit(1)
    
expr_dir=''
num_eq_classes=0
norm_file=''
unorm_file=''

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        expr_dir=arg
    elif opt in ("-m","--num-eq-classes"):
        num_eq_classes=int(arg)
    elif opt in ("-t","--TCC-file"):
        unorm_file=arg
    elif opt in ("-d","--TCCD-file"):
        norm_file=arg

if (not expr_dir) or (not num_eq_classes) or (not norm_file) or (not unorm_file):
    print ('usage is : \n python get_tcc_dist.py -i input_tcc_dir -m number-of-eq-classes -t path-to-output-TCC-file -d path-to-output-TCC-dist-file')
    sys.exit(1)

fl_list='file_list_Trapnell.dat'



# In[10]:


eq_dict={}

flnames=sorted([x for x in os.listdir(expr_dir) if x.endswith('.class')])
eq_class_hash=num_eq_classes

for flname in flnames:
    with open(expr_dir+flname) as flptr:
        for line in flptr:
            line = line.strip()
            vect = line.split()
            if not vect[0].isdigit():
                if vect[0] not in eq_dict:
                    eq_dict[vect[0]]=eq_class_hash
                    eq_class_hash+=1

TCC_mat=np.zeros((len(flnames),max(eq_dict.values())+1))


for cell_number in range(len(flnames)):
    cur_flname=flnames[cell_number]
    with open(expr_dir+cur_flname) as flptr1:
        for line in flptr1:
            line = line.strip()
            vect = line.split()
            assert len(vect)==2
            if vect[0].isdigit():
                index = int(vect[0])-1
            else:
                index = eq_dict[vect[0]]
            #print index
            value = int(vect[1])
            #print value
            TCC_mat[cell_number][index] = value

#print (np.shape(TCC_mat))
#print (sum(TCC_mat[0]>0))


TCC_dist, num_mapped_reads =data_to_dist(TCC_mat)
#print (TCC_dist.shape)
S=scipy.sparse.csr_matrix(TCC_dist)
S1=scipy.sparse.csr_matrix(TCC_mat)

with open(norm_file, 'wb') as outfile:
    pickle.dump(S, outfile, pickle.HIGHEST_PROTOCOL)

with open(unorm_file, 'wb') as outfile:
    pickle.dump(S1, outfile, pickle.HIGHEST_PROTOCOL)

with open(fl_list,'wb') as outfile:
    pickle.dump(flnames, outfile, pickle.HIGHEST_PROTOCOL)


# In[ ]:



