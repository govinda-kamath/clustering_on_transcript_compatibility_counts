from sklearn.metrics.pairwise import pairwise_distances
from scipy.stats import entropy
import pickle
import numpy as np
import sys
import multiprocessing as mp
import itertools

print len(sys.argv)
if len(sys.argv)!=4:
    print ('usage is \n python get_pairwise_distances.py ip-file op-file num-processes')
    exit(1)
#def jensen_shannon(pqtuple):
#    p=pqtuple[0]
#    q=pqtuple[1]
def jensen_shannon(p,q):
    #    pshape=np.shape(p)
    #    qshape=np.shape(q)
    #    assert pshape[1]==1
    #    assert qshape[1]==1
    #    assert pshape[0]==qshape[0]
    #    assert min(p) >= 0
    #    assert min(q) >= 0
    #assert sum(p)<=1 + np.finfo(float).eps and sum(p) >= 1- np.finfo(float).eps
    #assert sum(q)<=1 + np.finfo(float).eps and sum(q) >= 1- np.finfo(float).eps
    m=0.5*p+0.5*q
    p = np.transpose(p[p > 0])
    q = np.transpose(q[q > 0])
    m = np.transpose(m[m > 0])
     
    #if entropy(m)-0.5*entropy(q)-0.5*entropy(p) < 0:
    #    print (entropy(m)-0.5*entropy(q)-0.5*entropy(p))
    #    print ('m = '+ str(entropy(m)))
    #    print ('q = '+ str(entropy(q)))
    #    print ('p = '+ str(entropy(p)))
    return np.sqrt(entropy(m)-0.5*entropy(q)-0.5*entropy(p))

print (sys.argv[1])
print (sys.argv[2])
num_jobs=int(sys.argv[3])

with open(sys.argv[1], 'rb') as infile:
        Y = pickle.load(infile)
        print (np.shape(Y))

X=Y

num_rows=np.shape(X)[0]
for ind in range(num_rows):
    print(X[ind].sum(), ind)

D = pairwise_distances(X,metric=jensen_shannon,n_jobs=num_jobs)
#D=np.zeros((num_rows,num_rows))

#for i in range(num_rows):
#    print(i)
#    pool=mp.Pool(processes=num_jobs)
#    pqtuple=itertools.product([X[i,:]], X)
#    D[i,:]=pool.map(jensen_shannon,pqtuple)

with open(sys.argv[2],'wb') as outfile:
    pickle.dump(D, outfile, pickle.HIGHEST_PROTOCOL)


