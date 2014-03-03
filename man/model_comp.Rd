\name{model_comp}
\alias{model_comp}

\title{
Compare different models
}
\description{
This function makes it easy to compare different models based on different design matrices.
}
\usage{
model_comp(object, ...)
}
\arguments{
  \item{object}{
An object of class \code{nrm} or \code{nlm}, which is the outcome of fitting a nominal response model \code{nrm()} or a nested logit model \code{nelm()}.
}
  \item{\dots}{
Insert different design matrices here. Use the help function \code{designTemp()} to create different designs. These designs will be translated into an design matrix, and the supplied data will be fitted with these different designs.
}
}
\details{
The main advantage of this function is first, that different models with different design-matrices are fitted in one step and second that reasonable starting values are used, based on the values of the sumitted model object. This should accelerate the model fitting process (by reducing the number of iterations until convergence). \\

The loop which fits the models always \code{try()}s to fit the model, which means that an output is returned even if one or more models fail to converge.


}
\value{

A list of fitted models.

}

\author{
Manuel Reif
}



\seealso{
\code{\link{designTemp}}
\code{\link{reshMG}}
\code{\link{nrm}}
\code{\link{nelm}}

}
\examples{

\dontrun{
# Don't run this example because it takes more than 5sec.
NUMBI <- 5

ParList <- lapply(1:NUMBI,function(x)
{
  Item1 <- c(c(-2,-1,1,2),c(-1.2,0.3,0.2,0.7))
  names(Item1) <- c(paste("zeta",1:4,sep=""),paste("lamb",1:4,sep=""))
  Item1
})


names(ParList) <- paste("item",1:NUMBI,sep="")

set.seed(111)

perp1 <- rnorm(5000,0,1)
perp2 <- rnorm(5000,1,1)

simdat1 <- NRM.sim(ParList,perp1)
simdat2 <- NRM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:5000,simdat1)
simdat2 <- data.frame(ID=5001:10000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=TRUE)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=5000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)

reshOBJ <- reshMG(DAT1,items=1:NUMBI,groups=NUMBI+1,correct=rep(0,NUMBI))


# ---------------------------------------------------------------------------------

ergNRM <- nrm(reshOBJ)


mydes <- designTemp(ngru=2,nit=NUMBI,TYPE="NRM")
mydes[[1]][2,1] <- 2

erg_modc <- model_comp(ergNRM,mydes)


}

}

\keyword{model comparison}

