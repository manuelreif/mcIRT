\name{DDF}
\alias{DDF}

\title{
Differential Distractor Functioning
}
\description{
This function yields nonparametric measures for getting insight into items and response categories. There are two groups expected - a reference group and a focal group. The empirical item response curves are the basis of the resulting indices.
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
Which weights should be used? There are 3 options: "focal" (default), "reference" and "total". Focal is fine in the most cases. For further informations see the paper of Dorans and Kullick (1986).
}
}
\details{
\bold{reference group:} This denotes the group which is the "first" group within the reshape object - e.g. if the variable which denotes the different groups contains levels "A" and "B" then "A" is called the reference group and "B" is the focal group. You can choose the order of groups easily with \code{factor()} by permutating the levels. 

The choice which group to label as the reference group and which one as the focal group is up to the researcher. Often the larger group is chosen as the reference group, whereas the focal group is often the one which is of special interest. 


Only 2 groups are possible! If you want to compare > 2 groups, apply this function repeatedly.

}
\value{

\item{stdpdif}{Standardized difference for each category in each item}
\item{rmwsd}{RMWSD index}
\item{varpp}{Variance of the Probability Difference of reference minus focal group}
\item{Ps}{Points of the empirical item response function}

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

\code{\link{reshMG}}
\code{\link{lrt}}

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

# inspection of some arbitrary chosen item
my_ddf[[1]][[1]]


####### different reference group ##########

# use the level argument to set a different reference group
simdatallg2 <- data.frame(GROUP=factor(rep(c("A","B"),each=6000),levels=c("B","A")),simdatall)


##### RESHAPE ###########

my_resh2 <- reshMG(simdatallg2,items=2:11,groups=1,correct=rep(3,10),echo=FALSE)
  
##### DDF ###########

my_ddf2 <- DDF(my_resh2)

# inspection of some arbitrary chosen item
my_ddf2[[1]][[1]]


}



\keyword{DDF}

