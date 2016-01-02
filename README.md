# Clustering on transcript compatibility counts
Rather than clustering cells based on the abundances of genes, here we determine the transcript compatibility counts (TCCs) for each gene for each cell and cluster on the TCCs instead. We obtain transcript compatibility counts using [kallisto](https://github.com/pachterlab/kallisto).

To prove that using TCCs is a correct approach, we reanalyzed two recently published datasets:

1. The 271 primary human myoblasts by [Trapnell et al.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4122333/)
2. The 3005 mouse brain cells by [Zeisel et al.](http://linnarssonlab.org/cortex/)

See the corresponding folders for analysis code.
