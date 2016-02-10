import sys
import numpy as np
import os
import itertools
import pickle
import multiprocessing as mp
import timeit
import getopt

#This is created while running Zeisel_wrapper.py
file_list='./SRR_in_3005.txt'

cell_files=sorted(np.loadtxt(file_list,dtype=str))


try:
    opts, args = getopt.getopt(sys.argv[1:],"k:s:r:h:",["hacked-kallisto-path","seed=","reference-transcriptome","mouse-genome-for-hisat"])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python time_test.py -k hacked-kallisto-path -r mouse-reference-transcriptome -h path-to-file-containing-paths-to-mouse-genome [-s seed]')
    sys.exit(1)

kallipso_path=''
ref_transcriptome=''
ref_genome=''
npseed=100

for opt,arg in opts:
    if opt in ("-s", "--seed"):
        npseed=int(arg)
    elif opt in ("-k","--hacked-kallisto-path"):
        kallipso_path=arg
    elif opt in ("-r","--reference-transcriptome"):
        ref_transcriptome=arg
    elif opt in ("-h","--mouse-genome-for-hisat"):
        ref_genome=arg

if (not kallipso_path) or (not ref_transcriptome) or (not ref_genome):
    print ('usage is : \n python time_test.py -k hacked-kallisto-path -r mouse-reference-transcriptome -h path-to-file-containing-paths-to-mouse-genome [-s seed]')
    sys.exit(1)

test_SRA_dir='./SRA/'
test_read_dir='./reads/'
test_kallisto_dir='./kallisto/'
test_kallipso_dir='./TCC/'
test_bowtie1_dir='./bowtie1/'
test_hisat_dir='./hisat/'
test_wc_dir='./wc/'
kallisto_index='./kallisto_index/'
bowtie_index='./bowtie_index/'
hisat_index='./hisat_index/'

os.system('rm -rf '+test_SRA_dir)
os.system('rm -rf '+test_read_dir)
os.system('rm -rf '+test_kallisto_dir)
os.system('rm -rf '+test_kallipso_dir)
os.system('rm -rf '+test_bowtie1_dir)
os.system('rm -rf '+test_hisat_dir)
os.system('rm -rf '+test_wc_dir)
os.system('rm -rf '+kallisto_index)
os.system('rm -rf '+bowtie_index)
os.system('rm -rf '+hisat_index)

os.system('mkdir -p '+test_SRA_dir)
os.system('mkdir -p '+test_read_dir)
os.system('mkdir -p '+test_kallisto_dir)
os.system('mkdir -p '+test_kallipso_dir)
os.system('mkdir -p '+test_bowtie1_dir)
os.system('mkdir -p '+test_hisat_dir)
os.system('mkdir -p '+test_wc_dir)
os.system('mkdir -p '+kallisto_index)
os.system('mkdir -p '+bowtie_index)
os.system('mkdir -p '+hisat_index)

np.random.seed(seed=100)
files_picked=np.random.choice(cell_files,10,replace=False)


print('Copying over SRA files...')
base_cmd1="wget -O ./SRA/"
base_cmd2="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP045/SRP045452/"
for flnames in files_picked:
    cmd=base_cmd1+flnames+".sra "+base_cmd2+flnames+"/"+flnames+'.sra'
    os.system(cmd)
    
print('Converting SRA to fastq.gz files')
for flname in files_picked:
    cmd='fastq-dump --gzip '+ './SRA/'+flname+'.sra' ' -O '+ test_read_dir
    os.system(cmd)

print('Building kallisto index...')
kallisto_index_path=kallisto_index+'Zeisel_index.idx'
os.system(kallipso_path+' index -i '+kallisto_index_path+' '+ref_transcriptome)

print('Building hisat index...')
hisat_ip_paths=''
with open('hisat_chr_path_list.txt','r') as f:
    hisat_ip_paths=f.readline()
hisat_index_path=hisat_index+'Zeisel_index'
os.system('hisat-build  --offrate 5 '+hisat_ip_paths+' '+hisat_index_path)

print('Getting bowtie indices...')
bowtie_index_dir='./bowtie_index/'
os.system('mkdir -p '+bowtie_index_dir)
bowtie_index_path=bowtie_index_dir+'Zeisel_index.all'
os.system('bowtie-build --offrate=5 '+ref_transcriptome+' '+bowtie_index_path)


print('Timing kallisto...')
def run_kallisto():
    test_kallisto_dir='./kallisto/'
    test_read_dir='./reads/'
    transcriptome_path='./kallisto_index/Zeisel_index.idx'
    
    flnames=sorted(os.listdir(test_read_dir))
    for fls in flnames:
        cellname=fls.split('.')[0]
        out_dir=test_kallisto_dir+cellname
        read_fl=test_read_dir+fls
        cmd="""kallisto quant  -i"""+transcriptome_path+ """  -o """ +out_dir +""" --single -l 200 -s 100 """+read_fl
        #print cmd
        os.system(cmd)

x=timeit.timeit(run_kallisto,number=1)
op_file=test_kallisto_dir+'time.time'
with open(op_file,'w') as f:
    f.write(str(x))
    
print('Timing hacked kallisto...')
def make_kallipso_path_global():
    global kallipso_path
make_kallipso_path_global()
def run_kallipso():
    test_kallipso_dir='./TCC/'
    ref_path='./kallisto_index/Zeisel_index.idx'
    test_read_dir='./reads/'
    
    flnames=sorted(os.listdir(test_read_dir))
    for flname in flnames:
        read_path=test_read_dir+flname
        cellname=flname.split('.')[0]
        command = kallipso_path+' pseudoalign -i '+ ref_path+ ' -o ' + test_kallipso_dir+cellname+'.counts' + ' ' + read_path
        #print command
        os.system(command)
x=timeit.timeit(run_kallipso,number=1)
op_file=test_kallipso_dir+'time.time'
with open(op_file,'w') as f:
    f.write(str(x))
    
print('Timing word count...')
def run_word_count():
    test_read_dir='./reads/'
    flnames=sorted(os.listdir(test_read_dir))
    for flname in flnames:
        read_path=test_read_dir+flname
        os.system('wc '+read_path)
x=timeit.timeit(run_word_count,number=1)
op_file=test_wc_dir+'time.time'
with open(op_file,'w') as f:
    f.write(str(x))        

print('Timing hisat...')
def run_hisat():
    test_hisat_dir='./hisat/'
    test_read_dir='./reads/'
    transcriptome_path='./hisat_index/Zeisel_index'
    
    
    flnames=sorted(os.listdir(test_read_dir))
    for fls in flnames:
        cellname=fls.split('.')[0]
        out_sam=test_hisat_dir+cellname+'.sam'
        out_bam=test_hisat_dir+cellname+'.bam'
        read_fl=test_read_dir+fls
        cmd="""hisat -p 1 -x"""+transcriptome_path+ """  -U """ +test_read_dir+fls+" -S "+out_sam
        #print cmd
        cmd1="samtools view -bS "+out_sam+" > "+out_bam
        #print cmd1
        os.system(cmd)
        os.system(cmd1)
        
x=timeit.timeit(run_hisat,number=1)
op_file=test_hisat_dir+'time.time'
with open(op_file,'w') as f:
    f.write(str(x))

print('Timing bowtie...')
def run_bowtie1():
    index = './bowtie_index/Zeisel_index.all'
    test_bowtie1_dir='./bowtie1/'
    test_read_dir='./reads/'
    
    flnames=sorted(os.listdir(test_read_dir))
    for flname in flnames:
        out = test_bowtie1_dir + flname.split('.')[0]
        os.system('mkdir -p ' + out)
        BTcmd = 'gzip -dc '+test_read_dir+flname+' | bowtie  -p 1 -aS --offrate 1 '+index+' - | samtools view -Sb - > '+out+'/hits.bam'
        #print BTcmd
        os.system(BTcmd)
    
x=timeit.timeit(run_bowtie1,number=1)
op_file=test_bowtie1_dir+'time.time'
with open(op_file,'w') as f:
    f.write(str(x)) 

        