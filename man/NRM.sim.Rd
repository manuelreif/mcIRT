\name{NRM.sim}
\alias{NRM.sim}

\title{
Simulate data (nrm)
}
\description{
This function simulates data according to the Nominal Response Model.
}
\usage{
NRM.sim(Parlist, pp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Parlist}{A list of numerical vectors with item parameters has to be committed. Note that the zeta parameters as well as the lambda parameters have to add up to 0 for each item. For details concerning the structure of the list see section 'Details'.
}
  \item{pp}{
A numerical vector of person parameters.
}
}
\details{

The structure of the Parlist is as follows:

\itemize{

\item Vector for Item 1: c(zeta1, zeta2, \ldots, lambda1, lambda2, \ldots)
\item Vector for Item 2: c(zeta1, zeta2, \ldots, lambda1, lambda2, \ldots)
\item etc. 
}

}
\value{

The function returns a \code{length(pp)}*\code{length(Parlist)} data.frame.

}
\references{
Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored in two or
more nominal categories. \emph{Psychometrika, 37}, 29-51.
}
\author{
Manuel Reif
}


\seealso{
\code{\link{NLM.sim}}

\code{\link{nrm}}
}
\examples{

%\dontrun{

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

%}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nominal response model}
\keyword{simulate}
