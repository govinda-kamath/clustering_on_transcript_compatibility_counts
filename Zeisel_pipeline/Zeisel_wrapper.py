import os
import getopt
import sys


try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:k:t:",["idir=","njobs=","hacked-kallisto-path=","reference-transcriptome="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python Zeisel_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-human-reference-transcriptome [-n number-of-processes-to-use]')
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
    print ('usage is : \n python Zeisel_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-mouse-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)

print('Extracting reads from SRAs...')
os.system('mkdir -p ./reads_with_UMIs/')
os.system('rm -f ./reads_with_UMIs/*')
os.system('python process_SRA.py -i '+SRA_dir+' -o ./reads_with_UMIs/ -n '+str(num_proc))

print('Separating reads and UMIs...')
os.system('mkdir -p ./reads_and_UMI_subsample100/')
os.system('mkdir -p ./tmp_dir/')
os.system('rm -f ./reads_and_UMI_subsample100/*')
os.system('rm -f ./tmp_dir/*')
os.system('python Clean_reads.py -i ./reads_with_UMIs/ -o ./reads_and_UMI_subsample100/ '+
          '-t ./tmp_dir/ -n '+str(num_proc))
os.system('rmdir ./tmp_dir')

print('Generating the Kallisto index (with hacked kallisto)...')
os.system('mkdir -p ./kallisto_index')
os.system('rm -f ./kallisto_index/*')
index_path='./kallisto_index/Zeisel_index.idx'
os.system(kallipso_path+' index -i '+index_path+' '+ref_transcriptome)
metadata_cmd=kallipso_path+' metadata '+index_path
os.system(metadata_cmd)
num_ec = sum(1 for line in open('./kallisto_index/Zeisel_index.idx_ecmap.txt'))
print(num_ec)

print('Generating TCC (with hacked kallisto)...')
os.system('mkdir -p ./transcript_compatibility_counts_subsample100/')
os.system('rm -f ./transcript_compatibility_counts_subsample100/*')
os.system('python get_pseudoalignments.py -i ./reads_and_UMI_subsample100/ -o ./transcript_compatibility_counts_subsample100/ -k '+kallipso_path+ ' -t '+ index_path +' -n '+ str(num_proc))

print('Generating TCC distribution...')

os.system('python get_tcc_dist.py -i ./transcript_compatibility_counts_subsample100/ -m '+str(num_ec)+' -t ./Zeisel_TCC_subsample100.dat -d ./Zeisel_TCC_distribution_subsample100.dat')

print('Generating pairwise distances...')

os.system('python get_pairwise_distances.py ./Zeisel_TCC_distribution_subsample100.dat ./Zeisel_TCC_pairwise_JS_distance_subsample100.dat '+str(num_proc))


