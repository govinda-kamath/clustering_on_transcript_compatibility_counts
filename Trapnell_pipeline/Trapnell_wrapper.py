import os
import getopt
import sys


try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:",["idir=","njobs="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python Trapnell_wrapper.py -i input_SRA_dir [-n number-of-processes-to-use]')
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
    print ('usage is : \n python Trapnell_wrapper.py -i input_SRA_dir [-n number-of-processes-to-use]')
    sys.exit(1)

print('Extracting reads from SRAs...')
os.system('python process_SRA.py -i '+SRA_dir+' -o ./reads/ -n '+str(num_proc))