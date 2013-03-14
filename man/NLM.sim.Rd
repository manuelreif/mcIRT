\name{NLM.sim}
\alias{NLM.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates data according to the Nested Logit Model
}
\description{
This function simulates data according to the Nested Logit Model.
}
\usage{
NLM.sim(Parlist, pp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Parlist}{ A list of numerical vectors with item parameters has to be committed. Note that the zeta parameters as well as the lambda parameters have to add up to 0 for each item. For details concerning the structure of the list see section 'Details'.
}
  \item{pp}{
A numerical vector of person parameters.
}
}
\details{

The structure of the Parlist is as follows:

\itemize{

\item Vector for Item 1: c(slope, intercept, zeta1, zeta2, \ldots, lambda1, lambda2, \ldots)
\item Vector for Item 2: c(slope, intercept, zeta1, zeta2, \ldots, lambda1, lambda2, \ldots)
\item etc. 
}



}
\value{

The function returns a \code{length(pp)}*\code{length(Parlist)} data.frame.

}
\references{
Suh, Y., & Bolt, D. M. (2010). Nested logit models for multiple-choice item response data.
\emph{Psychometrika, 75}, 454-473.

}
\author{
Manuel Reif
}

\seealso{
\code{\link{NRM.sim}}

\code{\link{nelm}}
}
\examples{
\dontrun{
# item parameters

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


# person parameters
pp <- rnorm(1000)

# simulate
sim.nlm.1 <- NLM.sim(ParList,pp)
}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nested logit model}
\keyword{simulate}
