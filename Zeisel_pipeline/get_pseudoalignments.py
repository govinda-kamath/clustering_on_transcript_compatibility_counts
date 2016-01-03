
# coding: utf-8

# In[1]:

import os 
import itertools
import multiprocessing as mp
import getopt
import sys

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:o:n:k:t:",["idir=","odir=","njobs=","hacked-kallisto-path=","reference-transcriptome="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python get_pseudoalignments.py -i input-read-dir -o output-tcc-dir -k path-to-hacked-kallisto'+
           ' -t path-to-mouse-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)


read_dir=''
code_path=''
ref_path=''
out_dir=''
num_proc=1

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        read_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-k","--hacked-kallisto-path"):
        code_path=arg
    elif opt in ("-t","--reference-transcriptome"):
        ref_path=arg
    elif opt in ("-o","--odir"):
        out_dir=arg

# In[23]:
print ("Num processes = "+str(num_proc))

if (not read_dir) or (not code_path) or (not ref_path) or (not out_dir):
    print ('usage is : \n python get_pseudoalignments.py -i input-read-dir -o output-tcc-dir -k path-to-hacked-kallisto '+
           '-t path-to-mouse-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)

def get_pseudoalignments(fltuple):
    flname=fltuple[0]
    read_dir=fltuple[1]
    out_dir=fltuple[2]
    ref_path=fltuple[3]
    code_path=fltuple[4]
    out_path=out_dir+flname+'.class'
    read_path_1=read_dir+flname+'.fastq.gz'

    #print(flname)
    command = code_path+' pseudoalign -i '+ ref_path+ ' -o ' + out_path + ' ' + read_path_1
    #print(command)
    os.system(command)


# In[24]:

flnames=os.listdir(read_dir)
relevant_filenames=[files for files in flnames if files.endswith('.fastq.gz')]
#print relevant_filenames
file_names=sorted([x.split('.')[0] for x in relevant_filenames])
fltuples=itertools.product(file_names,[read_dir],[out_dir],[ref_path],[code_path])
#print file_names


# In[22]:

pool=mp.Pool(processes=num_proc)
#get_pseudoalignments(relevant_filenames[0])
pool.map(get_pseudoalignments,fltuples)

#297,305 Equivalence classes
