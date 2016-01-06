# Clustering on transcript compatibility counts
Rather than clustering cells based on their transcript abundances or gene expression, here we determine the transcript compatibility counts (TCCs) for each cell and cluster on the TCCs instead. We obtain transcript compatibility counts using [kallisto](https://github.com/pachterlab/kallisto).

To show that using TCCs is a correct approach, we reanalyzed two recently published datasets:

1. The 271 primary human myoblasts by [Trapnell et al.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4122333/)
2. The 3005 mouse brain cells by [Zeisel et al.](http://linnarssonlab.org/cortex/)

We obtained the raw read files from NCBI's [Gene Expression Omnibus](http://www.ncbi.nlm.nih.gov/geo/). See the corresponding folders for analysis code. The Trapnell_pipeline and Zeisel_pipeline folders contain scripts for automatically downloading the SRR files corresponding to datasets. We highly recommend looking at the documentation in the iPython notebooks Trapnell_Analysis.ipynb, Zeisel_Analysis.ipynb, and Timing_Analysis.ipynb. All notebooks go through how we use the code to generate the figures in our paper.

## Preliminaries

The following tools are required to run the pipelines in this repository.

* [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation)
* [Samtools](http://www.htslib.org/download/)
* [bowtie1](http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/)
* [bowtie2](http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/)
* [kallisto](https://github.com/pachterlab/kallisto)
* [eXpress](http://bio.math.berkeley.edu/eXpress/tutorial.html)
* [HISAT](https://ccb.jhu.edu/software/hisat/manual.shtml)

To run the scripts in the iPython notebooks, the following Python modules are required.

* [iPython](http://ipython.org/install.html)
* [scikit-learn](http://scikit-learn.org/stable/install.html)
* [numpy](http://docs.scipy.org/doc/numpy-1.10.1/user/install.html)
* [scipy](http://www.scipy.org/install.html)
* [matplotlib](http://matplotlib.org/users/installing.html#mac-osx-using-pip)
* [pysam](https://github.com/pysam-developers/pysam)
* [networkx](https://networkx.github.io/documentation/latest/install.html)

To run the code related to analysis on data of Trapnell et al., please follow the following instructions:

* Build the modified version of kallisto for paired end reads. This is in [modified-kallisto-paired](https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts/tree/master/modified-kallisto-source/kallisto_pseudo_paired). 
* Download the human transcriptome (for example from [here](http://bio.math.berkeley.edu/kallisto/transcriptomes/Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz)). 
* Download the data set of Trapnell et al. from [here](ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP033/SRP033135/) to get all the .sra files in a single directory. We've provided a sample script that can do this in [get_files.py](https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts/blob/master/Trapnell_pipeline/get_files.py). 
* Pass the directory of the SRA files, the path to the human trancriptome, and path to the modified version of kallisto for paired end reads to [Trapnell_Analysis.ipynb](https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts/blob/master/Trapnell_pipeline/Trapnell_Analysis.ipynb), to verify all results on the dataset of Trapnell et al.

To run the code related to analysis on data of Zeisel et al., please follow the following instructions:

* Build the modified version of kallisto for single ended reads. This is in [modified-kallisto-single](https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts/tree/master/modified-kallisto-source/kallisto_pseudo_single). This is also used for the timing analysis.
* Download the mouse transcriptome (for example from [here](http://bio.math.berkeley.edu/kallisto/transcriptomes/Mus_musculus.GRCm38.rel79.cdna.all.fa.gz))
* Download the data set of Trapnell et al. from [here](ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP045/SRP045452/) to get all the .sra files in a single directory. We've provided a sample script that can do this in [get_files.py](https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts/blob/master/Zeisel_pipeline/get_files.py). 
* Pass the directory of the SRA files, the path to the mouse trancriptome, and path to the modified version of kallisto for single end reads to [Zeisel_Analysis.ipynb](https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts/blob/master/Zeisel_pipeline/Zeisel_Analysis.ipynb), to verify all results on the dataset of Zeisel et al.


* Download the Mouse genome (for example from [here](ftp://ftp.ncbi.nih.gov/genomes/M_musculus/Assembled_chromosomes/seq/)) and gunzip all the fa.gz files. This is necessary only for running the timing analysis.

The figure below compares our TCC clustering pipeline to the conventional cell clustering pipeline. Please refer to our manuscript for more details.

![pipeline](https://github.com/govinda-kamath/clustering_on_reads/blob/master/pipeline.png)
