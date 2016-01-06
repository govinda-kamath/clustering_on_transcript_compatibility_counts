#!/usr/bin/env python

import os

os.system("mkdir -p ./SRA_files")
base_cmd1="wget -O ./SRA_files/SRR"
base_cmd2="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP045/SRP045452/SRR"
for index in xrange(1544693,1548184):
    cmd=base_cmd1+str(index)+".sra "+base_cmd2+str(index)+"/SRR"+str(index)+'.sra'
    os.system(cmd)
    

