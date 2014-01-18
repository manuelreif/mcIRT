\name{nelm}
\alias{nelm}
\alias{summary.nelm}
\alias{print.nelm}
\alias{deviance.nelm}

\title{
Estimating a Nested Logit Model
}
\description{

This function fits a Nested Logit model proposed by Suh and Bolt (2010). Like the Nominal Response Model this model is especially useful for multiple choice items. In contrast to the Nominal Response Model it models the correct answer category by means of a 2-PL model. For the distractors a NRM is fitted (for details take a closer look at the references mentioned below.)

}
\usage{
nelm(reshOBJ, etastart = "aut", ctrl=list())

\method{summary}{nelm}(object, \dots)
\method{print}{nelm}(x, \dots)
\method{deviance}{nelm}(object, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reshOBJ}{
An object of class reshNLM is expected. So the step before fitting the model is to reshape the data by means of the \code{reshMG} function.
}
  \item{etastart}{
A numerical vector or "aut" (which is the default). Starting values for the eta parameters can be changed (but is not necessary in typical cases).
}

  \item{ctrl}{
A list of argumets to customize the computations.
}

  \item{object}{
An object of class \code{nlm}.
}
  \item{x}{
An object of class \code{nlm}.
}

\item{\dots}{ 
}

}
\details{
The eta parameters in \code{etastart} denote the estimable parameters of the model. For example, for an item with 4 categories (1 correct answer and 3 distractors), 1 \eqn{\alpha}, 1 \eqn{\beta}, 2 \eqn{\gamma}'s (which substitute the 3 \eqn{\zeta}'s) and 2 \eqn{\xi}'s (which substitute the 3 \eqn{\lambda}'s) are constrained for the normalization (sum of parameter sets is zero). 



The following arguments can be comitted within a list (ctrl argument)

\itemize{
    \item \code{nodes}
  A numerical vector of length 1. Set the number of quadrature nodes for the a-priori distribution. The distribution is assumed to be normal.

    \item \code{absrange}
  A numerical vector of length 1. Denotes the absolute range of the a-priori distribution. The default value is 5, so the normal distribution ranges from \eqn{[-5 ; 5]}.
  
    \item \code{verbose}
  If TRUE, the estimation process is displayed in terms of the actual EMstep.
  
    \item \code{sigmaest}
  If TRUE, the variance of the latent person distribution is estimated. Otherwise it is set to 1 (for each group).
  
    \item \code{exac}
  A numerical vector of length 1. If the difference in the log-Likelihood between two consecutive EM steps is not larger than 'exac' - the estimation will stop. Default: \code{0.001}
  
    \item \code{EMmax}
  A numerical vector of length 1. This argument sets the maximum number of EM steps. The default value is 500. Feel free to enlarge this number.
  
  \item \code{NRmax}
A numerical vector of length 1. This argument sets the maximum number of Newton Raphson steps within the M-Step of the EM Algorithm. Default: \code{20}

\item \code{NRexac}
A numerical vector of length 1. If the difference between two consecutive NR steps is not larger than \code{NRexac} - the NR procedure stops. Default: \code{0.01}
  
\item \code{centBETA} If TRUE the estimated \eqn{\beta} parameters are centered to sum up to 0. Default is FALSE.
  
\item \code{centALPHA} If TRUE the estimated \eqn{\alpha} parameters are centered to sum up to 1. Default is FALSE.

\item \code{Clist}
A list which contains informations about which parameters should be held constant during estimation. Each list element has to look similar to this expression: "eta\\d* = -*\\d*" (of course real digits instead of Regexes!). The term \code{eta} refers to a column in the Q matrix in the reshape object which actually represents the parameters. The right side of the equation is the constant the parameter should be set to. So an entry in the \code{Clist} could look like "eta2 = -1" which means that the second eta parameter will not be estimated and is set to the value of -1.


\item \code{npar}
If TRUE, the prior distribution is nonparametric and is reestimated in each EM step by use of the expected number of examinees on each quadrature node. Default: FALSE. Computations are based on the EH (empirical histogram) estimation method of Woods 2007 & 2011. It is possible to estimate EHs for more than one group. First experiences showed that a huge amount of EM steps are needed (> 5000) to approximate the latent ability distributions.


\item \code{quads}
Supply specific quadrature nodes and weights for nonparametric estimation (\code{npar} must be TRUE). It has to be a list of length = number of groups. Inside each list element a list with two elements (a vector of nodes and a vector of weights) is expected. To get an idea what this looks like use \code{mcIRT:::quadIT(nodes=15,ngr=2)}.


  }

\value{


\item{etapar}{A numerical vector of eta-parameters}
\item{last_estep}{A list of informations concerning the last e-step before convergence. This is nothing the typical user should care about.}
\item{last_mstep}{The output provided by \code{optim} concerning the last M-step of the EM-Algorithm.}
\item{n_steps}{The number of passed EM steps.}
\item{erg_distr}{Estimates concerning the latent person distribution.}
\item{QUAD}{Denotes a list containing the quadrature nodes and weights which were used as a-priori distribution.}
\item{starting_values}{A list with infos concerning the starting values. The first entry gives merely the structure of the starting values whereas \code{$ulstv} gives the used starting values for the first EM step}
\item{EAPs}{The exact a-posteriori values for each person - which is a person parameter estimate. (Group membership is considered.)}
\item{ZLpar}{The list of item parameter estimates for each group.}
\item{SE}{The list of Standard Errors for the item parameter estimates.}
\item{reshOBJ}{The committed reshape object (which includes the data).}
\item{Catinf}{A list which contains 1) the information amount for each category/item/group for a sequence of ability values; 2) the sequence of ability values; 3) test information (sum above all items) for each of the ability values.}
\item{call}{Shows the actual call of the nrm function.}




}
\references{
Suh, Y., & Bolt, D. M. (2010). Nested logit models for multiple-choice item response data.
\emph{Psychometrika, 75}, 454-473.

Suh, Y. & Bolt, D. M. (2011). A nested logit approach for investigating distractors as causes of differential item functioning. \emph{Journal of Educational Measurement, 48}, 188-205.

Woods, C. M. (2007). Empirical Histograms in Item Response Theory With Ordinal Data. \emph{Education and Psychological Measurement}, 67:1, 73-87.

Woods, C. M. (2011). DIF Testing With an Empirical-Histogram Approximation of the Latent Density for Each Group. \emph{Applied Measurement in Education}, 24:3, 256-279.
}
\author{
Manuel Reif
}


\seealso{
\code{\link{reshMG}}

\code{\link{NLM.sim}}

\code{\link{nrm}}
}
\examples{

%\dontrun{
# create list of parameters
Item1 <- c(1,-2,c(-0.5,0.3,0.2),c(-0.5,-0.3,0.8))
names(Item1) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item2 <- c(1,-1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item2) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item3 <- c(1,0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item3) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item4 <- c(1,1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item4) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item5 <- c(1,2,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item5) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

ParList <- list(Item1=Item1,Item2=Item2,Item3=Item3,Item4=Item4,Item5=Item5)

# simulate data
perp1 <- rnorm(1000,0,1)
simdat1 <- NLM.sim(ParList,perp1)

# reshape
reshOBJ <- reshMG(simdat1,items=1:5,groups=NA,correct=rep(1,5),design="nodif",echo=TRUE,TYPE="NLM")

# estimate a nested logit model with a maximum number of 40 EM iterations, 
# which is NOT recommanded and is just applied here because estimating the model
# during example checks on cran took too long with default settings

res.nlm <- nelm(reshOBJ=reshOBJ, ctrl=list(EMmax=40))

summary(res.nlm)

%}

  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nested logit model}

