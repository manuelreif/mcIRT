\name{lrt}
\alias{lrt}
\alias{lrt.nrm}
\alias{lrt.nelm}
\alias{lrt.modc}


\title{
Likelihood Ratio Test
}
\description{
LRTs
}
\usage{
lrt(object, ...)
\method{lrt}{nrm}(object,  \dots)
\method{lrt}{nelm}(object,  \dots)
\method{lrt}{modc}(object,  \dots)

}

\arguments{
  \item{object}{
An object of class nrm, nelm or modc.
}
  \item{\dots}{
Designs
}
}
\details{
-- insert details here
}
\value{

Description of values.

}
\references{
REFS
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    UseMethod("lrt")
  }
}

\keyword{lrt}

