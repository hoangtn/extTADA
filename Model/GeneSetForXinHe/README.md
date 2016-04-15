We are using MCMC to estimate all parameters simultaneously, using this model:

$P(x| paramterers) = \prod_{i=1}^{m} \left[ \pi \textcolor{blue}{P(x_{i_{LoF}}|H_1)} \textcolor{brown}{P(x_{i_{mis3}}|H_1)} + (1 - \pi)\textcolor{blue}{P(x_{i_{LoF}}|H_0)} \textcolor{brown}{P(x_{i_{mis3}}|H_0)} \right] $


Missense damaging variants = missense variants + reported by 7 different methods

We are testing on these:

lof_maf001	AND	damaging_maf001

lof_Singleton	AND	missense_Singleton

lof_Singleton_noexac	AND	missense_Singleton_noexac

lof_maf001	AND	missense_maf001

lof_maf001_noexac	AND	missense_maf001_noexac

disruptive_Singleton_broad	AND	damaging_Singleton_broad

disruptive_Singleton_broad_noexac	AND	damaging_Singleton_broad_noexac

###########
### Some plots of current results are in *png files.

There are two pictures for Autism paper in 2014: 

1) Not use missense case-control in the estimation as the original paper: Tplot_lof_Singleton_noexac_missense_Singleton_noexac.png

2) Use missense case-control data in the estimation: Tplot_autism2014.png

############
### We intersected gene lists with other gene sets by:

1) Choose genes with FDR < 0.3 in:

intersect_with_giulio_genesets.ipynb

2) Choose top 100 genes from TADA results in:

intersect_with_giulio_genesets_top100.ipynb

Some significant results (overlapping genes) can be seen for FDR_scz2016_Singleon_noexac_pi026.txt.sort03 and FDR_scz2016_april.lof_maf001_missense_maf001pi.0.03.notMissenseInEstimation.txt 
