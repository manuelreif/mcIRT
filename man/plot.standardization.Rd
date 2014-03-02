\name{plot.standardization}
\alias{plot.standardization}

\title{
Plot DDF
}
\description{
This function visualizes DDF summary indices and some other more things in the future.
}
\usage{
\method{plot}{standardization}(x, type = "rmwsd", ...)
}

\arguments{
  \item{x}{
An object of class \code{standardization}.
}
  \item{type}{
Choose between different types of plots. Available: "rmwsd".
}
  \item{\dots}{
Further arguments for the \code{plot()} function.
}
}
\details{

There will be more plots ... soon.

}

\references{
Dorans, N. J. and Kullick, E. (1986). Demonstrating the utility of the standardization approach
to assessing unexpected differential item performance on the Scholastic Aptitude Test. Journal of 
Educational Measurement, 23, 355-368.

}
\author{
Manuel Reif
}



\seealso{
\code{\link{DDF}}
\code{\link{nrm}}
\code{\link{reshMG}}
}
\examples{

##### CREATE DATA ###########

x1 <- rnorm(10*5,0,1.7)
zeta <- tapply(x1,rep(1:10,each=5), function(x) round(x - mean(x),3))

# lambda parameters
x2 <- rnorm(10*5,0,1.1)
lam  <- tapply(x2,rep(1:10,each=5), function(x)
  {
  sort(round(x - mean(x),3),decreasing=FALSE)
  })


# create a parlist as a first step
ParList <- mapply(function(one,two)
  {
  x1 <- c(one,two)
  names(x1) <- paste0(rep(c("zeta","lam"),each=length(one)),1:length(one))
  return(x1)
  },one=zeta,two=lam,SIMPLIFY=FALSE)


names(ParList) <- paste0("item",1:length(ParList))


erglist     <- vector(mode="list",length=100)
erglistTIME <- vector(mode="list",length=100)

perp1 <- rnorm(6000,0,1)
perp2 <- rnorm(6000,0.3,1)

simdat1 <- NRM.sim(ParList,perp1)
simdat2 <- NRM.sim(ParList,perp2)

simdatall <- rbind(simdat1,simdat2)

simdatallg <- data.frame(GROUP=factor(rep(c("A","B"),each=6000)),simdatall)


##### RESHAPE ###########

my_resh <- reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE)
  
##### DDF ###########

my_ddf <- DDF(my_resh)


##### PLOT ###########

plot(my_ddf)


}

\keyword{DDF}

