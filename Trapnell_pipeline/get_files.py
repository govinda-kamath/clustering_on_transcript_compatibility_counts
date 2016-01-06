#!/usr/bin/env python

import os

os.system("mkdir -p ./SRA_files")
base_cmd1="wget -O ./SRA_files/SRR"
base_cmd2="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP033/SRP033135/SRR"
for index in xrange(1032910,1033294):
    cmd=base_cmd1+str(index)+".sra "+base_cmd2+str(index)+"/SRR"+str(index)+'.sra'
    os.system(cmd)

