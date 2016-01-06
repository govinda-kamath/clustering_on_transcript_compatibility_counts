
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
    return X


# In[3]:

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:t:d:",["idir=","TPM-file=","TPMD-file="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python get_UMI_count_matrices.py -i input_expression_dir  -d path-to-output-distribution-file -t path-to-unnormalised-count0-file')
    sys.exit(1)
    
expr_dir=''
norm_file=''
unorm_file=''

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        expr_dir=arg
    elif opt in ("-d","--TPMD-file"):
        norm_file=arg
    elif opt in ("-t","--TPM-file"):
        unorm_file=arg


if (not expr_dir) or  (not norm_file) or (not unorm_file):
    print ('usage is : \n python get_UMI_count_matrices.py -i input_expression_dir  -d path-to-output-distribution-file -t path-to-unnormalised-count0-file')
    sys.exit(1)

    
    


flnames=sorted([x for x in os.listdir(expr_dir) if x.endswith(".counts")])
test_fl_path=expr_dir+flnames[0]
num_transcripts = sum(1 for line in open(test_fl_path))

expr_mat=np.zeros((len(flnames),num_transcripts))

for ind in xrange(len(flnames)):
    flnm=expr_dir+flnames[ind]
    expr_mat[ind,:]=np.loadtxt(flnm, usecols=(1,))


normalised_expr_mat=data_to_dist(expr_mat)

S=scipy.sparse.csr_matrix(normalised_expr_mat)
S1=scipy.sparse.csr_matrix(expr_mat)

with open(norm_file, 'wb') as outfile:
    pickle.dump(S, outfile, pickle.HIGHEST_PROTOCOL)
with open(unorm_file, 'wb') as outfile:
    pickle.dump(S1, outfile, pickle.HIGHEST_PROTOCOL)