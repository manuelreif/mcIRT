\name{lrt}
\alias{lrt}
\alias{lrt.nrm}
\alias{lrt.nelm}
\alias{lrt.modc}


\title{
Likelihood Ratio Test
}
\description{
This function applies likelihood ratio tests, between one reference model (object) and some competing models.
}
\usage{
lrt(object, ...)
\method{lrt}{nrm}(object,  \dots)
\method{lrt}{nelm}(object,  \dots)
\method{lrt}{modc}(object,  \dots)

}

\arguments{
  \item{object}{
An object of class \code{nrm}, \code{nelm} or \code{modc}. In general, there are two approaches. The first one is fitting 2 or more models separately and supplying them to the function. The first object should be the most restrictive one (smalles number of parameters). The second approach is using the \code{model_comp()} function first to fit several models and then supplying the resulting object to this function to gain lrt tests.
}
  \item{\dots}{
Objects of class code{nrm} or \code{nelm}, which should be compared to the first object (reference model). So likelihood ratio tests are realized always between the reference model and all the others (1 vs 2, 1 vs 3, 1 vs 4 \ldots).
}
}
\details{
-- details? 
%% something about nesting structure when using lrt.
}
\value{

The function returns a matrix with the following values

\item{-2 loglikelihood difference}{This is the difference of every submitted model vs. the reference model (object).}
\item{degrees of freedom difference}{This is the difference in df between the reference model and the other models.}
\item{p-value}{The resulting X^2 value.}

}
\references{
Holland, P. W., & Wainer, H. (1993). Differential item functioning. Hillsdale, NJ: Lawrence Erlbaum.
}
\author{
Manuel Reif
}


\seealso{
\code{\link{nrm}}
\code{\link{nelm}}
\code{\link{model_comp}}
}
\examples{

## NUMBER OF ITEMS
nit <- 3

set.seed(2046)
x1 <- rnorm(nit*5,0,1.7)
zeta <- tapply(x1,rep(1:nit,each=5), function(x) round(x - mean(x),3))

# lambda parameters
x2 <- rnorm(nit*5,0,1.1)
lam  <- tapply(x2,rep(1:nit,each=5), function(x)
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


## NUMBER OF RESPONDENTS
nper <- 2000


names(ParList) <- paste0("item",1:length(ParList))

perp1 <- rnorm(nper,0,1)
perp2 <- rnorm(nper,0.3,1)

simdat1 <- NRM.sim(ParList,perp1)
simdat2 <- NRM.sim(ParList,perp2)

simdatall <- rbind(simdat1,simdat2)
simdatallg <- data.frame(GROUP=factor(rep(c("A","B"),each=nper)),simdatall)


my_resh <- reshMG(simdatallg,items=2:(nit+1),groups=1,correct=rep(3,nit),echo=FALSE)

est1 <- nrm(my_resh)


des2 <- designTemp(ngru=2,nit=nit,TYPE="NRM")
des2[[1]][2,1] <- 2 

my_resh2 <- reshMG(simdatallg,items=2:(nit+1),groups=1,correct=rep(3,nit),echo=FALSE,design=des2)

est2 <- nrm(my_resh2)

### LIKELIHOOD RATIO TEST
lrt(est1,est2)

}

\keyword{lrt}

