startV_nrmMG <-
function(reshOBJ,etastart=-0.1)
{

  granz <- length(levels(reshOBJ$gr))
  
  stwm <- lapply(1:granz,function(gy)
  {
    stwmi <- lapply(reshOBJ[[2]],function(x)
    {
      gam <- rep(0,length(x$categ))
      names(gam) <- rep("zet",length(x$categ))
      xi <- rep(0,length(x$categ))
      names(xi) <- rep("lam",length(x$categ))
      c(gam,xi)
    })
    stwmi
  })
  
  
  stwm1 <- as.relistable(stwm)          # skeleton
  ulstv <- vector(length=ncol(reshOBJ$Qmat),mode="numeric")
  ulstv[] <- etastart                   # values
  
  return(list(stwm1=stwm1,ulstv=ulstv))
}
