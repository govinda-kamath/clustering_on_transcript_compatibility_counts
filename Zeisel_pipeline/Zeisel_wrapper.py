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
#os.system('mkdir -p ./reads_with_UMIs/')
#os.system('rm -f ./reads_with_UMIs/*')
#os.system('python process_SRA.py -i '+SRA_dir+' -o ./reads_with_UMIs/ -n '+str(num_proc))

read_dir_base='./reads_and_UMI_subsample'
sampling_suffix=['100','10','5','1','_point5','_point1']

print('Separating reads and UMIs...')
read_dir_to_pass=read_dir_base+sampling_suffix[0]+"/"
#os.system('mkdir -p '+read_dir_to_pass)
#os.system('mkdir -p ./tmp_dir/')
#os.system('rm -f '+ read_dir_to_pass+'*')
#os.system('rm -f ./tmp_dir/*')
#os.system('python Clean_reads.py -i ./reads_with_UMIs/ -o '+read_dir_to_pass+' '+
#          '-t ./tmp_dir/ -n '+str(num_proc))
#os.system('rmdir ./tmp_dir')

print('Sampling reads')
sampling_rates=['1','0.1','0.05','0.01','0.005','0.001']
read_dir_to_pass=read_dir_base+sampling_suffix[0]+"/"
for index in range(1,6):
    print('Sampling '+sampling_rates[index]+' fraction of reads...')
    out_dir_to_pass=read_dir_base+sampling_suffix[index]+"/"
#    os.system('mkdir -p '+out_dir_to_pass)
#    os.system('rm -f '+out_dir_to_pass+'*')
    cmd='python sample_reads.py -i '+read_dir_to_pass+' -o '+out_dir_to_pass+' -k 100 -r '+sampling_rates[index]+' -n '+str(num_proc)
#    os.system(cmd)
    

print('Generating the Kallisto index (with hacked kallisto)...')
#os.system('mkdir -p ./kallisto_index')
#os.system('rm -f ./kallisto_index/*')
index_path='./kallisto_index/Zeisel_index.idx'
#os.system(kallipso_path+' index -i '+index_path+' '+ref_transcriptome)
#metadata_cmd=kallipso_path+' metadata '+index_path
#os.system(metadata_cmd)
num_ec = sum(1 for line in open('./kallisto_index/Zeisel_index.idx_ecmap.txt'))
print(num_ec)

print('Generating TCC (with hacked kallisto)...')
TCC_base_dir='./transcript_compatibility_counts_subsample'
for index in range(6):
    print('Running hacked kallisto on '+sampling_rates[index]+' fraction of reads...')
    TCC_dir=TCC_base_dir+sampling_suffix[index]+"/"
    read_dir_to_pass=read_dir_base+sampling_suffix[index]+"/"
#    os.system('mkdir -p '+TCC_dir)
#    os.system('rm -f '+TCC_dir+'*')
#    os.system('python get_pseudoalignments.py -i '+read_dir_to_pass+' -o '+TCC_dir+' -k '+kallipso_path+ ' -t '+ index_path +' -n '+ str(num_proc))

print('Generating TCC distribution...')
TCC_dist_base_flname='./Zeisel_TCC_distribution_subsample'
TCC_base_flname='./Zeisel_TCC_subsample'
for index in range(6):
    print('Getting the TCC dist for '+sampling_rates[index]+' fraction of reads...')
    TCC_dir=TCC_base_dir+sampling_suffix[index]+"/"
    TCC_dist_flname=TCC_dist_base_flname+sampling_suffix[index]+".dat"
    TCC_flname=TCC_base_flname+sampling_suffix[index]+".dat"
#    os.system('python get_tcc_dist.py -i '+TCC_dir+' -m '+str(num_ec)+' -t '+TCC_flname+' -d '+ TCC_dist_flname)

print('Generating pairwise distances...')
TCC_distance_base_flname='Zeisel_TCC_pairwise_JS_distance_subsample'
for index in range(6):
    TCC_dist_flname=TCC_dist_base_flname+sampling_suffix[index]+".dat"
    TCC_distance_flname=TCC_distance_base_flname+sampling_suffix[index]+".dat"
    print('Getting  pairwise distances for '+sampling_rates[index]+' fraction of reads...')
#    os.system('python get_pairwise_distances.py '+TCC_dist_flname+' '+TCC_distance_flname+' '+str(num_proc))

print('Running Kallisto quant and writing out pseudo-bams...')
quant_dir_base='Zeisel_kallisto_quant_subsample'
pbam_dir_base='Zeisel_pbam_subsample'
tmp_dir='./temp_dir/'
os.system('mkdir -p '+tmp_dir)
for index in range(6):
    print('Running kallisto pbam and quant for '+sampling_rates[index]+' fraction of reads...')
    pbam_dir= pbam_dir_base+sampling_suffix[index]+"/"
    quant_dir= quant_dir_base+sampling_suffix[index]+"/"
    read_dir_to_pass=read_dir_base+sampling_suffix[index]+"/"
    os.system('python run_kallisto_and_get_pbam.py -i '+read_dir_to_pass+' -o '+quant_dir+' -p '+pbam_dir+' -t '+tmp_dir+' -r '+index_path +' -n '+str(num_proc))
    
    
print('Getting kallisto matrices...')
kal_gene_dist_file_base='./Zeisel_kallisto_TPM_distribution_subsample'
for index in range(6):
    print('Getting kallisto matrices for '+sampling_rates[index]+' fraction of reads...')
    quant_dir= quant_dir_base+sampling_suffix[index]+"/"
    kal_gene_dist_file=kal_gene_dist_file_base+sampling_suffix[index]+".dat"
    os.system('python get_kallisto_matrices.py -i '+quant_dir+' -d '+kal_gene_dist_file)
