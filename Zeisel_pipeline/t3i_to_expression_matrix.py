# Obtain gene expression matrix for each of the samples sizes

import os
import sys
import numpy as np
import pickle
def data_to_dist(X):
    s = np.sum(X,axis=1)
    X = X / s[:,None]
    return X

def trans_exp_express(dirname):
    dict = {}
    files = np.array(sorted(os.listdir(dirname)))

    # First pass: get all unique transcripts, mapping each to an index
    count = 0
    for f in files:
        fo = open(dirname+'/'+f+'/results.t3i')
        for line in fo:
            temp = line.split()
            key = temp[0]
            if key not in dict: 
                dict[key] = count
                count += 1

    # Second pass: Count number of unique transcripts
    gene_abundances = np.zeros([len(files),len(dict)]);
    cell_ind = 0

    for f in files:
        fo = open(dirname+'/'+f+'/results.t3i')
        for line in fo:
            temp = line.split()
            key = temp[0]
            quant = float(temp[1])
            gene_abundances[cell_ind,dict[key]] += quant
        cell_ind += 1
        
    return gene_abundances

bowtie_base='./Zeisel_Bowtie_subsample'
norm_file_base='./Zeisel_express_distribution_subsample'
unorm_file_base='./Zeisel_express_subsample'
gene_file_base='./Zeisel_express_gene_subsample'
for suffix in ['100','10','5','1','_point5','_point1']:
    bowtie=bowtie_base+suffix
    unormX = trans_exp_express(bowtie)
    normX=data_to_dist(unormX)
    norm_file=norm_file_base+suffix+'.dat'
    gene_file=gene_file_base+suffix+'.dat'
    with open(norm_file,'wb') as outfile:
        pickle.dump(normX, outfile, pickle.HIGHEST_PROTOCOL)
    unorm_file=unorm_file_base+suffix+'.dat'
    with open(unorm_file,'wb') as outfile:
        pickle.dump(unormX, outfile, pickle.HIGHEST_PROTOCOL)

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

    gene_dist=np.zeros((np.shape(normX)[0],len(gene_hash)))
    col_ind=0
    for gid, cols in gene_hash.items():
        for colum in cols:
            gene_dist[:,col_ind]+=normX[:,colum]
        col_ind+=1



    with open(gene_file, 'wb') as outfile:
        pickle.dump(gene_dist, outfile, pickle.HIGHEST_PROTOCOL)

