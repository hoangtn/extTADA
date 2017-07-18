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


# Examples

Users can re-produce examples below to know how to use extTADA.

## Just one step to obtain final results

[./examples/extTADA_OneStep.ipynb]

## Some steps (to adjust plots for publication)

[./examples/extTADA_MultipleSteps.ipynb]


# Citation

If you use *extTADA*, please cite: 

*Bayesian Integrated Analysis Of Multiple Types Of Rare Variants To Infer Risk Genes For Schizophrenia And Other Neurodevelopmental Disorders*

Hoang T. Nguyen, Amanda Dobbyn, Laura M. Huckins, Douglas Ruderfer, Giulio Genovese, Menachem Fromer, Xinyi Xu, Joseph Buxbaum, Dalila Pinto, Christina Hultman, Pamela Sklar, Shaun M. Purcell, Xin He, Patrick F. Sullivan, Eli Ayumi Stahl

bioRxiv 135293; doi: https://doi.org/10.1101/135293

