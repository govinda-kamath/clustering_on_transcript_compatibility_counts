#! /bin/bash

(time python get_pairwise_distances_l1.py Zeisel_UMI_gene_distribution_subsample100.dat  Zeisel_UMI_test_l1.dat 32) > temp_op_UMI_l1.txt 2> time_UMI_l1.txt

(time python get_pairwise_distances_l1.py Zeisel_kallisto_TPM_distribution_subsample100.dat  Zeisel_kallisto_test_l1.dat 32) > temp_op_kall_l1.txt 2> time_kallisto_l1.txt

(time python get_pairwise_distances_l1.py Zeisel_TCC_distribution_subsample100.dat  Zeisel_TCC_test_l1.dat 32) > temp_op_TCC_l1.txt 2> time_TCC_l1.txt

(time python get_pairwise_distances.py Zeisel_UMI_gene_distribution_subsample100.dat  Zeisel_UMI_test.dat 32) > temp_op_UMI.txt 2> time_UMI.txt

(time python get_pairwise_distances.py Zeisel_kallisto_TPM_distribution_subsample100.dat  Zeisel_kallisto_test.dat 32) > temp_op_kall.txt 2> time_kallisto.txt

(time python get_pairwise_distances.py Zeisel_TCC_distribution_subsample100.dat  Zeisel_TCC_test.dat 32) > temp_op_TCC.txt 2> time_TCC.txt
