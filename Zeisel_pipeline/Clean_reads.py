# Separates the reads and UMIs in Zeisel's data. Takes care of the poly-G sequences.
# For more details, please refer to http://www.nature.com/nmeth/journal/v11/n2/full/nmeth.2772.html.

import itertools
import os
import multiprocessing as mp
import re
import getopt
import sys

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:o:n:t:",["idir=","odir=","njobs=","temp-dir="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python Clean_reads.py -i input-fastq-dir -o output-fastq-dir  -t temp-dir [-n number-of-processes-to-use]')
    sys.exit(1)
    
out_dir=''
read_dir=''
temp_dir=''
num_proc=1

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        read_dir=arg
    elif opt in ("-o","--odir"):
        out_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-t","--temp-dir"):
        temp_dir=arg
    
if (not read_dir) or (not temp_dir) or (not out_dir):
    print ('usage is : \n python Clean_reads.py -i input-fastq-dir -o output-fastq-dir  -t temp-dir [-n number-of-processes-to-use]')
    sys.exit(1)

def clean_file(fltuple):
    read_dir=fltuple[0]
    out_dir=fltuple[1]
    temp_dir=fltuple[2]
    flname=fltuple[3]
    cell_name=flname.split('.')[0]
    
    flpath=read_dir+flname
    cmd1= 'cp '+flpath+' '+temp_dir+flname
    cmd2= 'gunzip '+temp_dir+flname
    #cmd3= 'rm '+temp_dir+flname
    cmd5= 'rm '+temp_dir+cell_name+'.fastq'
    
    tmp_file=temp_dir+cell_name+'.fastq'
    out_file=out_dir+cell_name+'.fastq'
    umi_file=out_dir+cell_name+'.umi'
    
    cmd4='gzip '+out_file
    
    os.system(cmd1)
    os.system(cmd2)
    
    with open(tmp_file,'r') as f:
        with open(out_file,'w') as g:
            with open(umi_file,'w') as h:
                for line1,line2, line3, line4 in itertools.izip_longest(*[f]*4):

                    g.write(line1)
                    start = re.search(r'[^G]', line2[6:]).start()
                    g.write(line2[6+start:])
                    g.write(line3)
                    g.write(line4[6+start:])
                    h.write(line2[:6]+'\t' + str(start) +'\n')
    
    #os.system(cmd3)
    os.system(cmd4)
    os.system(cmd5)

flnames=sorted(os.listdir(read_dir))
fltuple=itertools.product([read_dir],[out_dir],[temp_dir],flnames)

pool=mp.Pool(processes=num_proc)
pool.map(clean_file,fltuple)
#for tup in fltuple:
#    clean_file(tup)

