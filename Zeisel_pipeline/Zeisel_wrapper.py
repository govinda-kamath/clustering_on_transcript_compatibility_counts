import os
import getopt
import sys


try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:",["idir=","njobs="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python Zeisel_wrapper.py -i input_SRA_dir [-n number-of-processes-to-use]')
    sys.exit(1)
    
SRA_dir=''
num_proc=1

for opt,arg in opts:
    #print (opt)
    #print (arg)
    if opt in ("-i", "--idir"):
        SRA_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
        
if not SRA_dir:
    print ('usage is : \n python Zeisel_wrapper.py -i input_SRA_dir [-n number-of-processes-to-use]')
    sys.exit(1)

print('Extracting reads from SRAs...')
#os.system('mkdir -p ./reads_with_UMIs/')
#os.system('rm -f ./reads_with_UMIs/*')
#os.system('python process_SRA.py -i '+SRA_dir+' -o ./reads_with_UMIs/ -n '+str(num_proc))

print('Separating reads and UMIs...')
os.system('mkdir -p ./reads_and_UMI_subsample100/')
os.system('mkdir -p ./tmp_dir/')
os.system('rm -f ./reads_and_UMI_subsample100/*')
os.system('rm -f ./tmp_dir/*')
os.system('python Clean_reads.py -i ./reads_with_UMIs/ -o ./reads_and_UMI_subsample100/ '+
          '-t ./tmp_dir/ -n '+str(num_proc))
os.system('rmdir ./tmp_dir')