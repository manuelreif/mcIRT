\name{plot.nrm}
\alias{plot.nrm}

\title{
Plotting an object of class nrm
}
\description{
This plotting routine draws the Category Characteristic Curves (CCC) for each item.
}
\usage{
\method{plot}{nrm}(x, numbpoints = 100, fromto = c(-4, 4), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Commit an object of class 'nrm' which is a result of the \code{nrm} function.}
  
  \item{numbpoints}{
A numerical vector of length 1, which denotes the number of points on the x-axis to draw the curves.
}
  \item{fromto}{
A numerical vector of length 2, which is about the same as xlim.
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
plot(res.nrm)

}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

