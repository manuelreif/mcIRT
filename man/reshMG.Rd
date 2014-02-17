\name{reshMG}
\alias{reshMG}

\title{
Reshaping the raw data
}
\description{
Run this function before fitting a nominal response model or nested logit model. The function prepares the data.frame for the estimation, adds and displays information and allows the user to submit a group estimation design which results in the creation of a Q-Matrix (Design Matrix).
}
\usage{
reshMG(da, items = 0, groups = NA, correct, design = "nodif",
echo = TRUE, TYPE = "NRM", paraM="bock")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{
The \code{data.frame} that has to be reshaped. Each row contains the responses of a person, each column contains the responses of different persons to one item. Use NA if a value is missing.
}
  \item{items}{
An index vector which indicates the columns containing the items. If it is set to 0 (which is the default), all columns are included.
}
  \item{groups}{
An index number which indicates the column containing the group membership of each person.
}
  \item{correct}{
An integer vector which denotes the position of the correct response category for each item. In cases in which no correct answer exists, it is reasonable to set this at \code{rep(0,length(items))}.
}
  \item{design}{
A character string or a group-design has to be committed. Valid inputs are: 'nodif' or a design as resulting from the \code{designTemp} function.
}
  \item{echo}{
Logical input, whether infos should be displayed on the console immediately after calling the function.
}
  \item{TYPE}{
Character string -- which kind of model the user intends to fit afterwards. There are two valid inputs: 'NRM' and 'NLM' 
}
  \item{paraM}{
Choose the parametrization. There are 2 valid inputs: "bock" and "01". If "bock" is used, the Q matrix is built following Bock 1972 - otherwise the last parameter of the NRM part is set so zero. To get a better idea of the differences take a closer look at the Q matrix after reshaping the data. 
}

}
\details{

The category responses in \code{da} must be \code{numeric} or \code{integer} and must start with 0. An item, for instance, with 5 options must be coded with 0,1,2,3,4. If e.g. the first category is the correct option set \code{correct} to 0 for this item.


In the case of a nominal response model, the \code{correct} argument just adds labels to the different categories, and does nothing else.

In the case of a nested logit model, the \code{correct} argument reorders the categories. The 'correct' option is set as first category and labeled as 'cor'. For example: If the vector c(1,0,2) is supplied, then for the first item the second (!) category is denoted as the 'correct' one. So inspecting the reshaped output we would see in the \code{recm} list element dummy matrices (in the case of more than 1 group), which are reordered due to the correct option (the correct option within each item is always first). The \code{d} element of the list, contains a list of recoded numeric matrices. The matrices were recoded so that '0' always denotes the 'correct' option. So in our example, the response vector for item 1 would be recoded so that each '1' would be recoded to '0' and each '0' to '1'. So \code{d} is in general not the same as the original input, but a version which is prepared for fitting the model.

The labels of the categories are consistent with the original numeric coding.


Note that a committed \code{design} is checked whether it is valid and possibly corrected. So it is recommended to check the function, the design and the design matrix.

}
\value{

 \item{recm}{A list with dummy matrices created from the (possibly transformed) original data.}
 \item{aDD}{Some infos about the data which are used by internal functions.}
 \item{d}{A list of response matrices - which may be different than the original input. See details for more information.}
 \item{gr}{A vector containing the group membership for each person.}
 \item{Qmat}{The design matrix. Feel free to modify the  desing matrix and insert it back into this object. Note that the \code{rownames} of the  design matrix \strong{must not be changed!}}
 \item{d1uc}{A modified version of the data input, which deals with NAs.}
 \item{design}{The committed design. The returned design might differ from the submitted one in order to get a valid design. Please check whether there were changes.}
}

\author{
Manuel Reif
}


\seealso{
\code{\link{designTemp}}

\code{\link{NRM.sim}}
}
\examples{

%\dontrun{

# create a list of parameters (see NRM.sim function)
NUMBI <- 5
ParList <- lapply(1:NUMBI,function(x)
{
  Item1 <- c(c(-2,-1,1,2),c(-1.2,0.3,0.2,0.7))
  names(Item1) <- c(paste("zeta",1:4,sep=""),paste("lamb",1:4,sep=""))
  Item1
})

names(ParList) <- paste("item",1:NUMBI,sep="")

# simulate the data for 2 groups
perp1 <- rnorm(1000,0,1)
perp2 <- rnorm(1000,1,1)

simdat1 <- NRM.sim(ParList,perp1)
simdat2 <- NRM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:1000,simdat1)
simdat2 <- data.frame(ID=1001:2000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=TRUE)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=1000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)

# reshape the data
reshdat <- reshMG(DAT1,items=1:5,groups=6,correct=rep(0,5), design="nodif")

# DIF design
mydes <- designTemp(ngru= 2, nit= 5, TYPE= "NRM")
mydes$zeta[2,5]   <- 2
mydes$lambda[2,5] <- 2

reshdat2 <- reshMG(DAT1,items=1:5,groups=6,correct=rep(0,5), design=mydes)

%}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{reshape}
\keyword{nominal response model}
\keyword{nested logit model}
