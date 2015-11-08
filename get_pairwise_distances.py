from sklearn.metrics.pairwise import pairwise_distances
from scipy.stats import entropy
import pickle
import numpy as np

def shannon_jensen(p, q):
    #    pshape=np.shape(p)
    #    qshape=np.shape(q)
    #    assert pshape[1]==1
    #    assert qshape[1]==1
    #    assert pshape[0]==qshape[0]
    #    assert min(p) >= 0
    #    assert min(q) >= 0
    #    assert sum(p)<=1 + np.finfo(float).eps and sum(p) >= 1- np.finfo(float).eps
    #    assert sum(q)<=1 + np.finfo(float).eps and sum(q) >= 1- np.finfo(float).eps
    m=0.5*p+0.5*q
    p = np.transpose(p[p > 0])
    q = np.transpose(q[q > 0])
    m = np.transpose(m[m > 0])
    return np.sqrt(entropy(m)-0.5*entropy(q)-0.5*entropy(p))

with open('/data/SS_RNA_seq/Code/sparse_eq_class.dat', 'rb') as infile:
        X = pickle.load(infile)
        print np.shape(X)

D = pairwise_distances(X,metric=shannon_jensen,n_jobs=64)

with open('/data/SS_RNA_seq/Code/pairwise_jensen_shannon_divergence.dat','wb') as outfile:
    pickle.dump(D, outfile, pickle.HIGHEST_PROTOCOL)


