# Clustering on transcript compatibility counts
Rather than clustering cells based on their transcript abundances or gene expression, here we determine the transcript compatibility counts (TCCs) for each cell and cluster on the TCCs instead. We obtain transcript compatibility counts using [kallisto](https://github.com/pachterlab/kallisto).

To show that using TCCs is a correct approach, we reanalyzed two recently published datasets:

1. The 271 primary human myoblasts by [Trapnell et al.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4122333/)
2. The 3005 mouse brain cells by [Zeisel et al.](http://linnarssonlab.org/cortex/)

See the corresponding folders for analysis code. We highly recommend looking at the documentation in the iPython notebooks Trapnell_Analysis.ipynb and Zeisel_Analysis.ipynb. Both notebooks go through how we use the code to generate the figures in our paper. 

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

The figure below compares our TCC clustering pipeline to the conventional cell clustering pipeline. Please refer to our manuscript for more details.

![pipeline](https://github.com/govinda-kamath/clustering_on_reads/blob/master/pipeline.png)
