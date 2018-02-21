# extTADA combines de novo mutations and case/control (transmitted/untransmitted) variants to:

1. Infer genetic parameters of these combined data sets.

2. Identify the number of significant genes based on the genetic parameters.

3. Predict number of significant genes based on a threshold using known genetic parameters (mainly from the tested data sets).

# extTADA can also be used:

1. For only de novo data OR only case/control (transmitted/untransmitted) data.

2. For multiple-population data.

## Requirement

These packages have to be installed:

1. *rstan*: [https://cran.r-project.org/web/packages/rstan/index.html](./https://cran.r-project.org/web/packages/rstan/index.html)

2. *locfit* [https://cran.rstudio.com/web/packages/locfit/index.html](./https://cran.rstudio.com/web/packages/locfit/index.html)

## Gene Lists

If you want to download top genes of schizophrenia (SCZ), intellectual disabilities (ID), developmental disorder (DD), epilepsy or autism spectrum disorder (ASD), you can use results analyzed by extTADA inside data. For example, top DD genes:

   *cat data/DDextTADA.txt |awk '$NF<0.1'*

OR download directly from the extTADA paper:

  * wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-017-0497-y/MediaObjects/13073_2017_497_MOESM2_ESM.xlsx*

# Examples 

Users can re-produce examples below to know how to use extTADA.

After running examples, please read [QandAfor_extTADA.md](./QandAfor_extTADA.md) if you have
questions.

## Just one step to obtain final results

[examples/extTADA_OneStep.ipynb](./examples/extTADA_OneStep.ipynb)

## Some steps (to adjust plots for publication)

[examples/extTADA_MultipleSteps.ipynb](./examples/extTADA_MultipleSteps.ipynb)

or (if cannot open ipython file)

[examples/extTADA_MultipleSteps.pdf](./examples/extTADA_MultipleSteps.pdf)


# Citation

If you use *extTADA*, please cite: 

*Integrated Bayesian analysis of rare exonic variants to identify risk genes for schizophrenia and neurodevelopmental disorders*

Hoang T. Nguyen, Julien Bryois, April Kim, Amanda Dobbyn, Laura M. Huckins, Ana B. Munoz-Manchado, Douglas M. Ruderfer, Giulio Genovese, 
Menachem Fromer, Xinyi Xu, Dalila Pinto, Sten Linnarsson, Matthijs Verhage, August B. Smit, Jens Hjerling-Leffler, Joseph D. Buxbaum, 
Christina Hultman, Pamela Sklar, Shaun M. Purcell, Kasper Lage, Xin He, Patrick F. Sullivan and Eli A. Stahl. Genome Medicine 2017 9:114.

Link: *https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0497-y#MOESM1*



