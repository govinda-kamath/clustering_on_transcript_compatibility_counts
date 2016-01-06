# coding: utf-8

# In[3]:

import sys
import pysam
import numpy as np
import os
import itertools
import pickle
import multiprocessing as mp
import getopt

# In[22]:
try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:k:t:",["idir=","njobs=","hacked-kallisto-path=","reference-transcriptome="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python UMI_counting.py  [-n number-of-processes-to-use]')
    sys.exit(1)
    
for opt,arg in opts:
    if opt in ("-n","--njobs"):
        num_proc=int(arg)

def count_UMI(celltuple):
    cellname=celltuple[0]
    cellpath=celltuple[1]
    print cellname
    transcript_path=celltuple[2]
    outfile=celltuple[3]
    umidir=celltuple[4]
    transcript_name=np.loadtxt(transcript_path,usecols=(0,),dtype=str)
    
    print "Reading bam file ..."
    bamfile_path=cellpath+cellname+'/'+'hits.bam'
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    umi_path=umidir+cellname+'.umi'
    
    read_start_alignment={}
    read_aligned_transcript={}
    read_id_set=set()
    mapped_reads=set()
    
    #print "Assigning reads to trancripts ..."
    for read in bamfile:
        read_id=int(read.query_name.split('.')[1])
        read_id_set.add(read_id)
        if not read.is_unmapped:
            mapped_reads.add(read_id)
            rd_start=int(read.get_reference_positions()[0])
            rd_trans=read.reference_name
            read_start_alignment.setdefault(read_id,rd_start)
            read_aligned_transcript.setdefault(read_id,rd_trans)
            if rd_start < read_start_alignment[read_id]:
                read_start_alignment[read_id]=rd_start
                read_aligned_transcript[read_id]=rd_trans
    
    #print "Reading UMIs and decoding their reads..."
    read_to_UMI={}
    UMI_list=np.loadtxt(umi_path,usecols=(0,),dtype=str)
    for read_id,index in itertools.izip(sorted(list(read_id_set)),xrange(len(UMI_list))):
        read_to_UMI[read_id]=UMI_list[index]
        
#     print len(read_id_set), len(UMI_list), len(mapped_reads)
#     return
    
    #print "Getting UMIs mapped to each transcript..."
    trans_to_umi={}
    transcripts_seen=read_aligned_transcript.values()
    for trans in transcripts_seen:
        trans_to_umi[trans]=set()
    for read_id in mapped_reads:
        trans_to_umi[read_aligned_transcript[read_id]].add(read_to_UMI[read_id])

    #print "Counting UMI on each transcript..."
    trans_to_num_UMI={}
    for trans in trans_to_umi:
        trans_to_num_UMI[trans]=len(trans_to_umi[trans])
    
    #print "Writing nonsense..."
    file_to_write=np.zeros(len(transcript_name),dtype=[('transname', 'S20'), ('UMIcount', 'S12')])
    for name,ind in itertools.izip(transcript_name, xrange(len(transcript_name))):
        file_to_write[ind]['transname']=name
        file_to_write[ind]['UMIcount']=str(trans_to_num_UMI.get(name, 0))
    
    np.savetxt(outfile+cellname+'.counts',file_to_write,fmt='%s')
    
    #print "DONE"


# In[23]:

cellpath_base='./Zeisel_Bowtie_subsample'
transcript_path='./transcript_length.txt'
outfile_base='./Zeisel_UMI_counts_subsample'
umi_dir='./reads_and_UMI_subsample100/'

sampling_suffix=['100']

for suffix in sampling_suffix:
    cellpath=cellpath_base+suffix+'/'
    outfile=outfile_base+suffix+'/'
    os.system('mkdir -p '+outfile)
    cellname=sorted([x for x in os.listdir(cellpath) if x.startswith('SRR')])
    celltuple=itertools.product(cellname,[cellpath],[transcript_path],[outfile],[umi_dir])

    pool=mp.Pool(processes=num_proc)
    pool.map(count_UMI,celltuple)