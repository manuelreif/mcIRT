\name{designTemp}
\alias{designTemp}

\title{
Create a design skeleton
}
\description{
This function creates a design skeleton -- which is only needed for multiple group data. This skeleton serves as template which can be changed by the user to define item by group interaction. The design is used as input for the reshape (\code{reshMG()}) function.
}
\usage{
designTemp(ngru, nit, TYPE = "NLM")
}

\arguments{
  \item{ngru}{
A digit which denotes the number of groups.
}
  \item{nit}{
A digit which denotes the number of items.
}
  \item{TYPE}{
Information about which model is going to be estimated. There are two valid values for this argument: 'NLM' and 'NRM' which refers to 'nested logit model' and 'nominal response model'.
}
}
\details{
For setting up a new multigroup design run \code{designTemp()} with default settings first. Subsequently modify the entries of the matrices in the template to get the user-defined multigroup-design.
The rows denote the different groups (e.g. interpreted as 'reference' and 'focal group'), and the columns denote the items (which depends on the chosen model). Each matrix in the list denotes a model specific parameter (zetas, lambdas etc.).


The first row in each matrix, has to contain only 1s, which means that these parameters will be estimated for the first group (e.g.: all the zeta parameters of item 1).


The second row in each matrix indicates if the parameter for this group is estimated seperately or together with another group. 
A '1' in the second row means that this parameter of group 2 is set equal to the parameter of the row above (which is in this example row 1 and therefore group 1).  
By contrast, a '2' in the second row means that no restriction is set to the parameter estimation. The parameter of group 2 is estimated independently from the parameter of group 1 
(because the first row is forced to contain only '1s'). For instance, if the zeta parameters of item 1 in group 2 should be estimated independently of those in group 1, fill in '2' (second row, first column) in the 'zeta' matrix. (To get an idea about the model parameters, take a closer look e.g. at Bock (1972), Baker & Kim (2004))


Perhaps there is a third group (therefore a third row in each matrix), and the zeta parameters of the first item should be constrained to be the same as in group 2, then fill in '2' in the third row/first column (zeta matrix). This means, that these parameters are estimated together with those of the group which is mentioned in the second row/ first column - which is, in our example, group 2.


}
\value{
The function returns a list of matrices.

}
\references{
Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored in two or
more nominal categories. \emph{Psychometrika, 37}, 29-51.

Baker, F. B., & Kim, S.-H. (2004). \emph{Item response theory. Parameter estimation techniques
(2nd ed.)}. New York, NY: Marcel-Dekker.
}
\author{
Manuel Reif
}


\seealso{
\code{\link{reshMG}}
}


\examples{

###### design 1 #####################

# DIF model - item 1 & 2 get their own zeta parameter estimates in group 2
des1 <- designTemp(2,5, "NRM")
des1$zeta[2,1:2] <- 2


###### design 2 #####################

# zetas: 
#item 5 - 9 are constrained to be the same in all 3 groups
#item 1: parameters are estimated in group 1 and 3 together, and separately for group 2
#item 2: same as item 1
#item 3: parameters are constrained to be the same in all 3 groups. Why? 
#Because the 2nd row/3rd column shows a '1' which means that the zeta parameters for the 
#second group are estimated together with the group mentioned in the 1st row/3rd column. 
#And 3rd row/ 3rd column shows a '2' which means, the third group is constrained to 
#estimate the same parameters as in the group which is mentioned in the 2nd row / 3rd column
#- which is the first group - so all parameters of all three groups are constrained to be the
#same for this item. Besides: it gives the same result as 3 '1' in a column.

#item 4: parameters are estimated in group 1 and group 2 together and separately in group 3
#lambdas:
# no DIF - all lambdas (for all items) are constrained to be the same in all groups.

des2 <- designTemp(3, 9, "NRM")
des2$zeta[2,1:2] <- 2
des2$zeta[3,3] <- 2
des2$zeta[3,4] <- 3


}
\keyword{nominal response model}
\keyword{nested logit model}