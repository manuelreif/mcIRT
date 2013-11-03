\name{plot_ldistr}
\alias{plot_ldistr}

\title{
Plot latent distributions
}
\description{
The resulting latent distribution is plotted. This is especially relevant when using a nonparametric estimation of the latent distribution.
}
\usage{
plot_ldistr(x,...)
}

\arguments{
  \item{x}{
An object of class \code{nrm} or \code{nelm}.
}

  \item{...}{
Insert further arguments to the plot function.
}
}
\details{

The function is plotting quadrature nodes and weights for each group (in the case of a multigroup estimation) which are approximating the posterior ability distribution for the group.

}

\author{
Manuel Reif
}


\seealso{
\code{\link{nrm}}

\code{\link{nelm}}

\code{\link{plotINF}}

}
\examples{

# Simulating 3 Items within a loop
ParList <- lapply(1:3,function(x)
{
  Item1 <- c(c(-2,-1,1,2),c(-1.2,0.3,0.2,0.7))
  names(Item1) <- c(paste("zeta",1:4,sep=""),paste("lamb",1:4,sep=""))
  Item1
})

names(ParList) <- paste("item",1:3,sep="")

# simulate person parameters
perp1     <- rnorm(1500,0,1)

# simulate data.frame 
sim.nrm.1 <- NRM.sim(ParList,perp1)

# reshape
reshdat1 <- reshMG(sim.nrm.1,items=1:3,correct=rep(1,3))

# fit the nrm
res.nrm <- nrm(reshdat1)

# finally the plot
plot_ldistr(res.nrm)



}

\keyword{plot}
