import sys
import pysam
import numpy as np
import os
import itertools
import pickle
import multiprocessing as mp
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:],"s:n:b:o:",["subsampled-read-dir=","njobs=","base-bam-dir=","out-bam-dir="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python subsample_bams.py -s subsampled-read-dir -b base-bam-dir -o out-bam-dir [-n number-of-processes-to-use]')
    sys.exit(1)
    

num_proc=1
bamfile_base_dir=''
output_base_dir=''
sub_reads_base_dir = ''

for opt,arg in opts:
    if opt in ("-s", "--subsampled-read-dir"):
        sub_reads_base_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-b","--base-bam-dir"):
        bamfile_base_dir=arg
    elif opt in ("-o","--out-bam-dir"):
        output_base_dir=arg

if (not bamfile_base_dir) or (not output_base_dir) or (not sub_reads_base_dir):
    print ('usage is : \n python subsample_bams.py -s subsampled-read-dir -b base-bam-dir -o out-bam-dir [-n number-of-processes-to-use]')
    sys.exit(1)
# In[5]:

def subsample_bam(celltuple):
    cellname=celltuple[0]
    bamfile_base_dir=celltuple[1]
    output_base_dir=celltuple[2]
    sub_reads_base_dir = celltuple[3]
    bamfile_path=bamfile_base_dir+cellname+'/'+'hits.bam'
    
    #bamfile_path='/data/SS_RNA_seq/Zeisel/samtools_experiments/little.bam'
    print cellname
#    print '\n'+'loading bamfile...'
    ## INPUT 1 -- original bamfile
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
#    print 'sorting bam...'
    read_dict={}
    for read in bamfile:
        if not read.is_unmapped:
            read_dict.setdefault(int(read.query_name.split('.')[1]),[])
            read_dict[int(read.query_name.split('.')[1])].append(read)
        
    sampling_suffix=['10','5','1','_point5', '_point1']

    
    for suffix in sampling_suffix:
#        print 'in '+suffix+'...'
        
        out_dir=output_base_dir+suffix+'/'
        os.system('mkdir -p '+ out_dir)
        out_cell_dir=out_dir+cellname+'/'
        os.system('mkdir -p '+ out_cell_dir)
        outfile_path=out_cell_dir+'hits.bam'
        
        sub_read_dir=sub_reads_base_dir+suffix+'/'
        sub_reads_path=sub_read_dir+cellname+'.rnum'
        umi_path=sub_read_dir+cellname+'.umi'
        sub_reads = sorted(np.loadtxt(sub_reads_path,dtype=int))
        umi=np.loadtxt(umi_path,dtype=str,usecols=(0,))
        
#        print 'writing header...'
        ## OUTPUT -- new bam file with the same header
        outfile = pysam.AlignmentFile(outfile_path, "wb", header=bamfile.header)
        
        mapped_UMI=[]
        unmapped_UMI=[]
#        print 'subsampling reads...'
        for index in xrange(len(sub_reads)):
            read_id=sub_reads[index]
            try:
                for read in read_dict[read_id]: 
                    outfile.write(read)
                mapped_UMI.append(umi[index])
            except KeyError:
                unmapped_UMI.append(umi[index])
                
        outfile.close()
        
#        print 'writing UMIs ...'
        
        np.savetxt(out_cell_dir+'mapped.umi',np.array(mapped_UMI),delimiter="\n", fmt="%s")
        np.savetxt(out_cell_dir+'unmapped.umi',np.array(unmapped_UMI),delimiter="\n", fmt="%s")
        
#         with open(out_cell_dir+'mapped.umi', 'wb') as f:
#             pickle.dump(mapped_UMI,f,pickle.HIGHEST_PROTOCOL)
#         with open(out_cell_dir+'unmapped.umi', 'wb') as f:
#             pickle.dump(unmapped_UMI,f,pickle.HIGHEST_PROTOCOL)
#        print 'DONE.'


# In[6]:

cellnames=sorted([x for x in os.listdir(bamfile_base_dir) if x.startswith('SRR')])

celltuples=itertools.product(cellnames,[bamfile_base_dir],[output_base_dir],[sub_reads_base_dir])
pool=mp.Pool(processes=num_proc)
pool.map(subsample_bam,celltuples)
