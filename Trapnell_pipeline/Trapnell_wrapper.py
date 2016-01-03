import os
import getopt
import sys
import numpy as np

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:k:t:",["idir=","njobs=","hacked-kallisto-path=","reference-transcriptome="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python Trapnell_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-human-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)
    
SRA_dir=''
num_proc=1
kallipso_path=''
ref_transcriptome=''

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        SRA_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-k","--hacked-kallisto-path"):
        kallipso_path=arg
    elif opt in ("-t","--reference-transcriptome"):
        ref_transcriptome=arg
    
        
if (not SRA_dir) or (not kallipso_path) or (not ref_transcriptome):
    print ('usage is : \n python Trapnell_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-human-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)

print('Extracting reads from SRAs...')
os.system('mkdir -p ./reads/')
os.system('rm -f ./reads/*')
os.system('python process_SRA.py -i '+SRA_dir+' -o ./reads/ -n '+str(num_proc))

print('Removing files that Trapnell throws...')
relevant_fls=np.loadtxt('Files_to_keep.txt',dtype=str)
fls_to_remove=[cells for cells in os.listdir('./reads/') if cells.split('_')[0] not in relevant_fls ]
for fl in fls_to_remove:
    os.system('rm '+'./reads/'+fl)


print('Generating the Kallisto index (with hacked kallisto)...')
os.system('mkdir -p ./kallisto_index')
os.system('rm -f ./kallisto_index/*')
index_path='./kallisto_index/Trapnell_index.idx'
os.system(kallipso_path+' index -i '+index_path+' '+ref_transcriptome)
metadata_cmd=kallipso_path+' metadata '+index_path
os.system(metadata_cmd)
num_ec = sum(1 for line in open('./kallisto_index/Trapnell_index.idx_ecmap.txt'))
print(num_ec)

print('Generating TCC (with hacked kallisto)...')
os.system('mkdir -p ./transcript_compatibility_counts/')
os.system('rm -f ./transcript_compatibility_counts/*')
os.system('python get_pseudoalignments_paired_end.py -i ./reads/ -o ./transcript_compatibility_counts/ -k '+kallipso_path+ ' -t '+ index_path +' -n '+ str(num_proc))

print('Generating TCC distribution...')

os.system('python get_tcc_dist.py -i ./transcript_compatibility_counts/ -m '+str(num_ec)+' -t ./Trapnell_TCC.dat -d ./Trapnell_TCC_distribution.dat')

print('Generating pairwise distances...')

os.system('python get_pairwise_distances.py ./Trapnell_TCC_dist.dat ./Trapnell_TCC_pairwise_distance.dat '+str(num_proc))
