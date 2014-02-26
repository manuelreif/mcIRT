\name{mcIRT-package}
\alias{mcIRT-package}
\alias{mcIRT}
\docType{package}
\title{
IRT models for multiple choice items (mcIRT)
}
\description{
This package provides functions to evaluate multiples choice items or other multicategorial response formats. It is possible to estimate multiple group models and to model interaction effects to examine Differential Item Functioning. Each model has a full accessible design matrix which allows the user to manipulate and set up his own weighting scheme and his own constraints. One application could be modeling and testing item properties by means of customizing the design matrix and estimating an explanatory Nominal Response Model or an explanatory Nested Logit Model (to get an idea what 'explanatory' means, take a look at Boeck & Wilson (2004)).
}
\details{
\tabular{ll}{
Package: \tab mcIRT\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2014-02-27\cr
License: GPL (>= 2)
}

}
\author{
Manuel Reif

Maintainer: Manuel Reif <manuel.reif@univie.ac.at>

}
\references{
Baker, F. B., & Kim, S.-H. (2004). \emph{Item response theory. Parameter estimation techniques
(2nd ed.)}. New York, NY: Marcel-Dekker.

Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored in two or more nominal categories. \emph{Psychometrika, 37}, 29-51.

Suh, Y., & Bolt, D. M. (2010). Nested logit models for multiple-choice item response data.
Psychometrika, 75, 454-473.

Suh, Y. & Bolt, D. M. (2011). A nested logit approach for investigating distractors as causes of differential item functioning. \emph{Journal of Educational Measurement, 48}, 188-205. 

De Boeck, P. & Wilson, M. (2004). \emph{Explanatory Items Response Models}. New York, NY:
Springer.

Woods, C. M. (2007). Empirical Histograms in Item Response Theory With Ordinal Data. \emph{Education and Psychological Measurement}, 67:1, 73-87.

Woods, C. M. (2011). DIF Testing With an Empirical-Histogram Approximation of the Latent Density for Each Group. \emph{Applied Measurement in Education}, 24:3, 256-279.

}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
%~~ directory ~~
\keyword{package}

\seealso{
\code{\link{nrm}}

\code{\link{nelm}}
}
\examples{

\dontrun{

# Simulating 5 Items within a loop
ParList <- lapply(1:5,function(x)
{
  Item1 <- c(c(-2,-1,1,2),c(-1.2,0.3,0.2,0.7))
  names(Item1) <- c(paste("zeta",1:4,sep=""),paste("lamb",1:4,sep=""))
  Item1
})

names(ParList) <- paste("item",1:5,sep="")

# simulate person parameters
perp1     <- rnorm(1500,0,1)

# simulate data.frame 
sim.nrm.1 <- NRM.sim(ParList,perp1)

# reshape
reshdat1 <- reshMG(sim.nrm.1,items=1:5,correct=rep(0,5))

# fit the nrm
res.nrm <- nrm(reshdat1)

}


}
