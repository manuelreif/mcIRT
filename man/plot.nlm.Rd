\name{plot.nelm}
\alias{plot.nelm}

\title{
Plotting an object of class nlm
}
\description{
This plotting routine draws the Category Characteristic Curves (CCC) for each item.
}
\usage{
\method{plot}{nelm}(x, numbpoints = 100, fromto = c(-4, 4), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Commit an object of class 'nlm' which is a result of the \code{nelm} function.
}
  \item{numbpoints}{
A numerical vector of length 1, which denotes the number of points on the x-axis to draw the curves.
}
  \item{fromto}{
A numerical vector of length 2, which is about the same as xlim.
}
  \item{\dots}{
More plot arguments.
}
}

\author{
Manuel Reif
}

\seealso{
\code{\link{nelm}}
}
\examples{
\dontrun{

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

# simulate
perp1 <- rnorm(1000,0,1)
simdat1 <- NLM.sim(ParList,perp1)

# reshape
reshOBJ <- reshMG(simdat1,items=1:5,groups=NA,correct=rep(1,5),design="nodif",echo=TRUE,TYPE="NLM")

# estimate a nested logit model
res.nlm <- nelm(reshOBJ=reshOBJ)

# plot the estimated data
plot(res.nlm)


}
}
\keyword{plot}
