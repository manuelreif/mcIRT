\name{plotINF}
\alias{plotINF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Information Curves
}

\description{
Plot Category/Item/Test Information functions for nominal response models and nested logit models.
}
\usage{
plotINF(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{nrm} or \code{nelm}.
}

  \item{\dots}{
More arguments for plot.
}
}

\author{
Manuel Reif
}


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
reshdat1 <- reshMG(sim.nrm.1,items=1:5,correct=rep(1,5))

# fit the nrm
res.nrm <- nrm(reshdat1)

# finally the plot
plotINF(res.nrm)
}
}
\keyword{plot}
