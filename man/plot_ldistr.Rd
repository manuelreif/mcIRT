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

\dontrun{
## simulating data

NUMBI <- 5
ParList <- lapply(1:NUMBI,function(x)
{
  Item1 <- c(c(-2,-1,1,2),c(-1.2,0.3,0.2,0.7))
  names(Item1) <- c(paste("zeta",1:4,sep=""),paste("lamb",1:4,sep=""))
  Item1
})


names(ParList) <- paste("item",1:NUMBI,sep="")

perp1 <- rnorm(5000,0,1)
perp2 <- rnorm(5000,1,1)

simdat1 <- NRM.sim(ParList,perp1)
simdat2 <- NRM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:5000,simdat1)
simdat2 <- data.frame(ID=5001:10000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=T)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=5000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)

## reshaping the data
reshOBJ <- reshMG(DAT1,items=1:NUMBI,groups=NUMBI+1,correct=rep(1,NUMBI))


## estimating the model
examp1 <- nrm(reshOBJ, ctrl=list(nodes=31))

} 

# loading the estimated model, because estimating the model during example check takes too long

data(examp1)

plot_ldistr(examp1)



}

\keyword{plot}
