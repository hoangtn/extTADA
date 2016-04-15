We are using MCMC to estimate all parameters simultaneously, using this model:


$P(x| paramterers) = \prod_{i=1}^{m} \left[ \pi \textcolor{blue}{P(x_{i_{LoF}}|H_1)} \textcolor{brown}{P(x_{i_{mis3}}|H_1)} + (1 - \pi)\textcolor{blue}{P(x_{i_{LoF}}|H_0)} \textcolor{brown}{P(x_{i_{mis3}}|H_0)} \right] $


Missense damaging variants = missense variants + reported by 7 different methods

We are testing on these:

lof_maf001	damaging_maf001
lof_Singleton	missense_Singleton
lof_Singleton_noexac	missense_Singleton_noexac
lof_maf001	missense_maf001
lof_maf001_noexac	missense_maf001_noexac
disruptive_Singleton_broad	damaging_Singleton_broad
disruptive_Singleton_broad_noexac	damaging_Singleton_broad_noexac


We intersected gene lists with other gene sets by:

1) Choose genes with FDR < 0.3 in:

intersect_with_giulio_genesets.ipynb

2) Choose top 100 genes from TADA results in:

intersect_with_giulio_genesets_top100.ipynb
