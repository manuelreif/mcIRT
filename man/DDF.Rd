\name{DDF}
\alias{DDF}

\title{
Differential Distractor Functioning
}
\description{
This function yields nonparametric measures for getting insight into items and distractors. There are two groups expected - a reference group and a focal group. The empirical item response curves are the basis of the resulting indices.
}
\usage{
DDF(reshOBJ, wm = "focal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reshOBJ}{
A reshape object which is a result of the \code{reshMG()} function. Make sure there are two groups supplied - the first is always supposed to be the reference group.
}
  \item{wm}{
Which weights should be used? There are 3 options: "focal" (default), "reference" and "total". Focal is fine in the most cases. For further informations see the references.
}
}
\details{

No more detail until now.

}
\value{

\item{stdpdif}{Standardized difference for each category in each item}
\item{rmwsd}{RMWSD index}
\item{varpp}{Variance of the Probability Difference of reference minus focal group}
\item{Ps}{Points of the empirical item response function}

}
\references{
More references
}
\author{
Manuel Reif
}


\seealso{

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

}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DDF}

