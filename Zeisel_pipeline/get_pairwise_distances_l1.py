from sklearn.metrics.pairwise import pairwise_distances
from scipy.stats import entropy
import pickle
import numpy as np
import sys
import multiprocessing as mp
import itertools

print(len(sys.argv))
if len(sys.argv)!=4:
    print ('usage is \n python get_pairwise_distances.py ip-file op-file num-processes')
    exit(1)
#def jensen_shannon(pqtuple):
#    p=pqtuple[0]
#    q=pqtuple[1]
print ("l1 distance")
print (sys.argv[1])
print (sys.argv[2])
num_jobs=int(sys.argv[3])

with open(sys.argv[1], 'rb') as infile:
        X = pickle.load(infile)
        print (np.shape(X))

#Y=X

def gk_manhattan(p,q):
    intermed=p-q
    sgn=intermed.sign()
    return intermed.dot(sgn.T).sum()

#num_rows=np.shape(Y)[0]
#for ind in range(num_rows):
#    print(X[ind].sum(), ind)

D = pairwise_distances(X,metric=gk_manhattan,n_jobs=num_jobs)
#D=np.zeros((num_rows,num_rows))

#for i in range(num_rows):
#    print(i)
#    pool=mp.Pool(processes=40)
#    pqtuple=itertools.product([X[i,:]], X)
#    D[i,:]=pool.map(jensen_shannon,pqtuple)

with open(sys.argv[2],'wb') as outfile:
    pickle.dump(D, outfile, pickle.HIGHEST_PROTOCOL)


