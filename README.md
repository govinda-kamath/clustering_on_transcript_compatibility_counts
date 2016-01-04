# Clustering on transcript compatibility counts
Rather than clustering cells based on their transcript abundances or gene expression, here we determine the transcript compatibility counts (TCCs) for each cell and cluster on the TCCs instead. We obtain transcript compatibility counts using [kallisto](https://github.com/pachterlab/kallisto).

To show that using TCCs is a correct approach, we reanalyzed two recently published datasets:

1. The 271 primary human myoblasts by [Trapnell et al.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4122333/)
2. The 3005 mouse brain cells by [Zeisel et al.](http://linnarssonlab.org/cortex/)

See the corresponding folders for analysis code. We highly recommend looking at the documentation in the iPython notebooks Trapnell_Analysis.ipynb and Zeisel_Analysis.ipynb. Both notebooks go through how we use the code to generate the figures in our paper. 

The figure below compares our TCC clustering pipeline to the conventional cell clustering pipeline. Please refer to our manuscript for more details.

![pipeline](https://github.com/govinda-kamath/clustering_on_reads/blob/master/pipeline.png)
