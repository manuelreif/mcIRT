\name{nrm}
\alias{nrm}
\alias{summary.nrm}
\alias{print.nrm}
\alias{deviance.nrm}

\title{
Estimating a Nominal Response Model
}
\description{
This function fits a Nominal Response Model as proposed by Bock (1972).
The model estimates two parameters for each response category: \eqn{\zeta} and \eqn{\lambda}. The model ist best suited for a set of items with more than two nominal response categories.

\deqn{P_{iqv} = \frac{exp(\zeta_{ik} - \lambda_{ik} \theta_j)}{\sum_{v=1}^{m} exp(\zeta_{iv} - \lambda_{iv} \theta_j)}}


}
\usage{
nrm(reshOBJ,etastart=-0.1, ctrl=list())
\method{summary}{nrm}(object, \dots)
\method{print}{nrm}(x, \dots)
\method{deviance}{nrm}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reshOBJ}{
A object of class reshNRM is expected. So the step before fitting the model is to reshape the data by means of the \code{reshMG} function.
}
  \item{etastart}{
A numerical vector. Starting values for the eta parameters can be change (but is not necessary in typical cases).
}
  \item{ctrl}{
A list of argumets to customize the computations.
}

  \item{object}{
An object of class \code{nrm}.
}

  \item{x}{
An object of class \code{nrm}.
}
\item{\dots}{
}
}
\details{

The eta parameters in \code{etastart} denote the estimable parameters of the model. For example, for an item with 4 categories, 3 \eqn{\gamma}'s (which substitute the 4 \eqn{\zeta}'s) and 3 \eqn{\xi}'s (which substitute the 4 \eqn{\lambda}'s) are constrained for normalization (parameter sets add up to zero). 


The following arguments can be comitted within a list (ctrl argument):

\itemize{
  \item \code{nodes}
A numerical vector of length 1. Set the number of quadrature nodes/points for the a-priori distribution. The distribution is assumed to be normal. Default: \code{14}

  \item \code{absrange}
A numerical vector of length 1. Denotes the absolute range of the a-priori distribution. The default value is 5, so the normal distribution ranges from \eqn{[-5 ; 5]}.

  \item \code{verbose}
If TRUE, the estimation process is displayed in terms of the actual EMstep. Default: \code{TRUE}

  \item \code{sigmaest}
If TRUE, the variance of the latent person distribution is estimated. Otherwise it is set to 1 (for each group). Default: \code{FALSE}

  \item \code{exac}
A numerical vector of length 1. If the difference between two consecutive EM steps (sum of absolute difference of parameters) is not larger than 'exac' - the estimation stops. Default: \code{0.00001}

  \item \code{EMmax}
A numerical vector of length 1. This argument sets the maximum number of EM steps. The default value is 500. Feel free to enlarge this number. Default: \code{500}

\item \code{NRmax}
A numerical vector of length 1. This argument sets the maximum number of Newton Raphson steps within the M-Step of the EM Algorithm. Default: \code{20}

\item \code{NRexac}
A numerical vector of length 1. If the difference between two consecutive NR steps is not larger than 'exac' - the estimation stops. Default: \code{0.01}

}


}
\value{

\item{etapar}{A numerical vector of eta-parameters}
\item{last_estep}{A list of informations concerning the last e-step before convergence. This is nothing the typical user should care about.}
\item{last_mstep}{Contains the log likelihood and the Hessian-Matrix of the last M-step of the EM-Algorithm.}
\item{n_steps}{The number of passed EM steps which lead to the convergent result.}
\item{erg_distr}{Estimates concerning the latent person distribution.}
\item{QUAD}{Denotes a list of the quadrature nodes and weights which were used as a-priori distribution.}
\item{starting_values}{A list with infos concerning the starting values. The first entry gives merely the structure of the starting values whereas \code{$ulstv} gives the used starting values for the first EM step}
\item{EAPs}{The exact a-posteriori values for each person - which is something like a person parameter estimate.}
\item{ZLpar}{The list of item parameter estimates for each group.}
\item{SE}{The list of standard errors for the item parameter estimates.}
\item{reshOBJ}{The committed reshape object (which includes the data).}
\item{Catinf}{A list which contains 1) the information amount for each category/item/group for a sequence of ability values; 2) the sequence of ability values; 3) test information (sum above all items) for each of the ability values.}
\item{call}{Shows the actual call of the nrm function.}


}
\references{
Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored in two or
more nominal categories. \emph{Psychometrika, 37}, 29-51.

Baker, F. B., & Kim, S.-H. (2004). \emph{Item response theory. Parameter estimation techniques
(2nd ed.)}. New York, NY: Marcel-Dekker.

Harwell, M. R., Baker, F. B. & Zwarts, M. (1988). Item parameter estimation via marginal maximum likelihood and an EM algorithm: A didactic. \emph{Journal of Educational Statistics, 13}, 243-271.

}
\author{
Manuel Reif
}


\seealso{
\code{\link{reshMG}}

\code{\link{NRM.sim}}

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
reshdat1 <- reshMG(sim.nrm.1,items=1:5,correct=rep(1,5))

# fit the nrm
res.nrm <- nrm(reshdat1)

summary(res.nrm)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nominal response model}
