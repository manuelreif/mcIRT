mcIRT
=====

This is an R package to fit 2 Item Response Theory (IRT) models: the Nominal Response Model and the Nested Logit Model.


**Version 0.40 will be released in the next few days on cran.**
This version is already available here!

To install this package, install [devtools](https://github.com/hadley/devtools) first.

```R
library(devtools)
install_github("mcIRT", "manuelreif", ref="master")
```


TO-DOs for the next release version:


* add some new examples and extend existing examples [creating examples in wiki]
* add info about convergence [done]
* fit a bunch of models in a row with smartly chosen starting values [done]
* do some simulations [done - see wiki - upcoming more]
* Nonparametric tests [first function with plot]
* add an LRT routine which computes Likelihood ratio tests for two or more models [done]
* watch the number of EM-steps growing until convergence ;-) [done]


Further TO-DOs:


* perhaps: add the possibility to resume estimation after `EMmax` is reached
* Big data adaption
* more tests :-)


If there are some features you want to see in future versions, don`t hesitate to write an email or report on github.


