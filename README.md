mcIRT
=====

This is an R package to fit 2 Item Response Theory (IRT) models, namely the Nominal Response Model and the Nested Logit Model.


To install this package, install devtools first.

```R
library(devtools)
install_github("mcIRT", "manuelreif", ref="master")
```

To install the development version, which now includes c++ code (thanks to the `Rcpp` package developers) for time consuming procedures install the development branch **nrmcpp**:

```R
library(devtools)
install_github("mcIRT", "manuelreif", ref="nrmcpp")
```

This version runs stable in first tests and provides markedly faster estimation of the nominal response model and the nested logit model as well (since 18.01.2014).

TO-DOs for the next release version:


* add some new examples and extend existing examples [creating examples in wiki]
* add info about convergence [done]
* add an LRT routine which computes Likelihood ratio tests for two or more models
* perhaps: add the possibility to resume estimation after `EMmax` is reached
* fit a bunch of models in a row with smartly chosen starting values [done]
* do some simulations [done - see wiki - upcoming more in the next days]
* watch the number of EM-steps growing until convergence [done]


Further TO-DOs:


* Nonparametric tests
* Big data adaption
* more tests :-)


If there are some features you want to see in future versions, don`t hesitate to write an email or report on github.


