
# coding: utf-8

# In[35]:

import os 
import subprocess
import itertools
import multiprocessing as mp
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:o:n:p:t:r:",["idir=","odir=","njobs=","pbam-path=","temp-dir=","ref-index-path="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python run_kallisto_and_get_pbam.py -i input-read-dir -o output-kallisto-dir -p pbam-dir'+
           ' -t path-to-temp-dir -r path-to-mouse-reference-index  [-n number-of-processes-to-use]')
    sys.exit(1)


out_dir=''
out_pbam_dir=''
read_dir=''
temp_dir=''
ref_path=''
num_proc=1

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        read_dir=arg
    elif opt in ("-o","--odir"):
        out_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-p","--pbam-path"):
        out_pbam_dir=arg
    elif opt in ("-r","--ref-index-path"):
        ref_path=arg
    elif opt in ("-t","--temp-dir"):
        temp_dir=arg

        
if (not out_dir) or (not out_pbam_dir) or (not read_dir) or (not temp_dir) or (not ref_path):
    print ('usage is : \n python run_kallisto_and_get_pbam.py --bias -i input-read-dir -o output-kallisto-dir -p pbam-dir'+
           ' -t path-to-temp-dir -r path-to-mouse-reference-index  [-n number-of-processes-to-use]')
    sys.exit(1)

# In[36]:

def get_pseudobam(fltuple):
    flname=fltuple[0]
    read_dir=fltuple[1]
    ref_path=fltuple[2]
    out_dir_base=fltuple[3]
    out_bam_dir=fltuple[4]
    temp_dir=fltuple[5]
    
    cellname=flname.split('.')[0]
    outdir=out_dir_base+cellname+'/'
    read_fl=read_dir+flname
    temp_file=temp_dir+cellname+'.bam'
    outfl=out_bam_dir+cellname+'.pbam'
    
    cmd="""kallisto quant -i """+ref_path+""" -o """ +outdir +""" --bias --single -l 200 -s 100 --pseudobam """+read_fl+" > "+temp_file
    cmd2="""samtools view -S """+temp_file  +""" | awk '{split($1,a,"."); split($12,b,":"); print a[2], $3, $4,b[3]}' > """+outfl
    cmd3='rm -f '+temp_file

    print cmd
    print cmd2
    print cmd3


    os.system(cmd)
    os.system(cmd2)
    os.system(cmd3)


# In[40]:



os.system('mkdir -p '+out_dir)
os.system('mkdir -p '+out_pbam_dir)

flnames=sorted([x for x in os.listdir(read_dir) if x.endswith('.fastq.gz')])
fltuple=itertools.product(flnames,[read_dir],[ref_path],[out_dir],[out_pbam_dir],[temp_dir])

pool=mp.Pool(processes=num_proc)
pool.map(get_pseudobam,fltuple)
#for tup in fltuple:
#    get_pseudobam(tup)

