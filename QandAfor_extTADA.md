# Questions 


## How can I input for multiple populations + categories?

You should use populations as categories. For example, if
you have two case-control populations and two categories then you
should use as four categories.

Example: 

First popultation: Ncase1 = 2000, Ncontrol1 = 3000

Second population: Ncase2 = 1000, Ncontrol2 = 5000

Then Ncase should be: Ncase = c(Ncase1, Ncase1, Ncase2, Ncase2)

and Ncontrol should be: Ncontrol = c(Ncontrol1, Ncontrol1, Ncontrol2,
Ncontrol2)

And the input data will have: 

cc_case1 = case counts from the *first* category of the *first* population

cc_case2 = case counts from the *second* category of the *first* population

cc_case3 = case counts from the *first* category of the *second* population

cc_case4 = case counts from the *second* category of the *second* population

A similar way should be used for control data.

## Can we manually use other distributions and/or numerical integration to calculate Bayes Factors (e.g,., not negative binomial distributions)?

Yes, people can do that easily. We have a note here:

https://github.com/hoangtn/information_from_other_papers/

We also tested results from the "dnbinom" function (used inside the original TADA model) and from numerical integration here:
https://github.com/hoangtn/information_from_other_papers/blob/master/file/Compare_BayesFactors.ipynb





