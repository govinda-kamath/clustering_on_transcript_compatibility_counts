import sys
import numpy as np
import os
import itertools
import pickle
import multiprocessing as mp
import timeit
import getopt

#This is created while running 
full_Zeisel_read_dir='../Zeisel_pipeline/reads_and_UMI_subsample100/'

read_files=sorted([x for x in os.listdir(full_Zeisel_read_dir) if x.endswith('.fastq.gz')])
files_picked=np.random.choice(read_files,10,replace=False)
files_used_in_paper=['SRR1545085.fastq.gz', 'SRR1546711.fastq.gz', 'SRR1547170.fastq.gz',
                     'SRR1547184.fastq.gz',  'SRR1547881.fastq.gz', 'SRR1545232.fastq.gz',
                     'SRR1547119.fastq.gz',  'SRR1547175.fastq.gz', 'SRR1547531.fastq.gz',  
                     'SRR1548035.fastq.gz']

try:
    opts, args = getopt.getopt(sys.argv[1:],"k:pr:h:",["hacked-kallisto-path","as-in-paper=","reference-transcriptome","mouse-genome-for-hisat"])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python time_test.py -k hacked-kallisto-path -r mouse-reference-transcriptome -h path-to-file-containing-paths-to-mouse-genome [-p] \n -p indicates to use the same files as used in the simulation of the paper')
    sys.exit(1)

kallipso_path=''
ref_transcriptome=''
ref_genome=''

for opt,arg in opts:
    if opt in ("-p", "--as-in-paper"):
        files_picked=files_used_in_paper
    elif opt in ("-k","--hacked-kallisto-path"):
        kallipso_path=arg
    elif opt in ("-r","--reference-transcriptome"):
        ref_transcriptome=arg
    elif opt in ("-h","--mouse-genome-for-hisat"):
        ref_genome=arg

if (not kallipso_path) or (not ref_transcriptome) or (not ref_genome):
    print ('usage is : \n python time_test.py -k hacked-kallisto-path -r mouse-reference-transcriptome -h path-to-file-containing-paths-to-mouse-genome [-p] \n -p indicates to use the same files as used in the simulation of the paper')
    sys.exit(1)
    
test_read_dir='./reads/'
test_kallisto_dir='./kallisto/'
test_kallipso_dir='./TCC/'
test_bowtie1_dir='./bowtie1/'
test_hisat_dir='./hisat/'
test_wc_dir='./wc/'
kallisto_index='./kallisto_index/'
bowtie_index='./bowtie_index/'
hisat_index='./hisat_index/'

os.system('rm -rf '+test_read_dir)
os.system('rm -rf '+test_kallisto_dir)
os.system('rm -rf '+test_kallipso_dir)
os.system('rm -rf '+test_bowtie1_dir)
os.system('rm -rf '+test_hisat_dir)
os.system('rm -rf '+test_wc_dir)
#os.system('rm -rf '+kallisto_index)
#os.system('rm -rf '+bowtie_index)
#os.system('rm -rf '+hisat_index)
os.system('mkdir -p '+test_read_dir)
os.system('mkdir -p '+test_kallisto_dir)
os.system('mkdir -p '+test_kallipso_dir)
os.system('mkdir -p '+test_bowtie1_dir)
os.system('mkdir -p '+test_hisat_dir)
os.system('mkdir -p '+test_wc_dir)
os.system('mkdir -p '+kallisto_index)
os.system('mkdir -p '+bowtie_index)
os.system('mkdir -p '+hisat_index)


print('Building kallisto index...')
kallisto_index_path=kallisto_index+'Zeisel_index.idx'
#os.system(kallipso_path+' index -i '+kallisto_index_path+' '+ref_transcriptome)

print('Building hisat index...')
hisat_ip_paths=''
with open('hisat_chr_path_list.txt','r') as f:
    hisat_ip_paths=f.readline()
hisat_index_path=hisat_index+'Zeisel_index'
os.system('hisat-build  --offrate 5 '+hisat_ip_paths+' '+hisat_index_path)