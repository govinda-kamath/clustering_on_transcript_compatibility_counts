
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
    opts, args = getopt.getopt(sys.argv[1:],"i:t:d:",["idir=","num-eq-classes=","TPM-file=","TPMD-file="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python get_kallisto_matrices.py -i input_expression_dir  -d path-to-output-distribution-file')
    sys.exit(1)
    
expr_dir=''
norm_file=''


for opt,arg in opts:
    if opt in ("-i", "--idir"):
        expr_dir=arg
    elif opt in ("-d","--TPMD-file"):
        norm_file=arg

if (not expr_dir) or  (not norm_file) :
    print ('usage is : \n python get_kallisto_matrices.py -i input_expression_dir  -d path-to-output-distribution-file')
    sys.exit(1)

    
    


flnames=sorted([x for x in os.listdir(expr_dir) if x.startswith("SRR")])
test_fl_path=expr_dir+flnames[0]+'/abundance.tsv'
num_transcripts = sum(1 for line in open(test_fl_path))-1

expr_mat=np.zeros((len(flnames),num_transcripts))

for ind in xrange(len(flnames)):
    flnm=expr_dir+flnames[ind]+'/abundance.tsv'
    expr_mat[ind,:]=np.loadtxt(flnm, skiprows=1, usecols=(4,))


normalised_expr_mat=data_to_dist(expr_mat)

trans_list=np.loadtxt('./transcript_length.txt',usecols=(0,),dtype=str) 
trans_hash={}
for trans_id, trans in enumerate(trans_list):
    trans_hash[trans]=trans_id

gene_hash={}
with open ('./transcript_to_gene_map.txt','r') as f:
    for line in f:
        trans_id=line.split()[0]
        gene_id=line.split()[1]
        gene_hash.setdefault(gene_id,set())
        gene_hash[gene_id].add(trans_hash[trans_id])

gene_dist=np.zeros((np.shape(normalised_expr_mat)[0],len(gene_hash)))
col_ind=0
for gid, cols in gene_hash.items():
    for colum in cols:
        gene_dist[:,col_ind]+=normalised_expr_mat[:,colum]
    col_ind+=1

S=scipy.sparse.csr_matrix(gene_dist)


with open(norm_file, 'wb') as outfile:
    pickle.dump(S, outfile, pickle.HIGHEST_PROTOCOL)




# In[ ]:

