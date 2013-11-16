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

This version runs stable in first tests and provides markedly faster estimation of the nominal response model. The c++ conversion of the nested logit model will follow soon!

