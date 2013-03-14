Cnlm <-
function(reshOBJ,ESTlist)
{
  
  retrans1 <- as.vector(reshOBJ$Qmat %*% ESTlist$etapar)
  
  catzjegr <- sapply(reshOBJ$aDD,function(x)x$anz_cat)
  catanz <- rep(rep(catzjegr,each=2),each=nlevels(reshOBJ$gr))-1
  
  # beta - parametrisierung 1
  wobeta <- grep("beta",rownames(reshOBJ$Qmat))
  supv1 <- rep(levels(reshOBJ$gr),each=length(reshOBJ$aDD))
  ergbet <- tapply(retrans1[wobeta],supv1,function(x) x -mean(x))
  
  woalpha <- grep("alpha",rownames(reshOBJ$Qmat))
  #supv1 <- rep(levels(reshOBJ$gr),each=length(reshOBJ$aDD))
  ergal <- tapply(retrans1[woalpha],supv1,function(x) x + (1 - mean(x)))
  
  # beta - parametrisierung 2
  beta2   <- retrans1[wobeta] / retrans1[woalpha]
  ergbet2 <- tapply(beta2,supv1,function(x) x -mean(x))
  
  # create mats
  matList <- mapply(function(scr,zae)
  {

      toD <- matrix( - 1/scr,ncol=scr,nrow=scr)
      diag(toD) <- diag(toD) + 1

    toD
  },scr=catanz,zae=1:length(catanz),SIMPLIFY=FALSE)  
  
  
  
  to2 <- cumsum(catanz)
  from1   <- c(1,to2[-length(to2)]+1)
  
  cem1 <- matrix(0,ncol=sum(catanz),nrow=sum(catanz))
  
  
  for(i in 1:length(from1))
  {
    fr <- from1[i]
    tt <- to2[i]
    
    cem1[fr:tt,fr:tt] <- matList[[i]]
  }
  
  nrmpartres <- as.vector(retrans1[-c(wobeta,woalpha)] %*% cem1)
  names(nrmpartres) <- rownames(reshOBJ$Qmat)[-c(wobeta,woalpha)]
  
  # dummylist
  sch2 <- lapply(levels(reshOBJ$gr),function(nooneknows)
            {
            sch1 <- lapply(reshOBJ$aDD,function(runfun)
                    {
                    zetas   <- rep(0,(runfun$anz_cat - 1))
                    lambdas <- rep(0,(runfun$anz_cat - 1))
                    list(zetas=zetas,lambdas=lambdas)
                    })
            names(sch1) <- paste("Item",1:length(sch1),sep="")
            sch1
            })
  names(sch2) <- levels(reshOBJ$gr)
  
  ergrell <- relist(nrmpartres,sch2)
  
  
  list("alpha"=ergal,"beta"=ergbet,"beta_differentPar"=ergbet2,"nrmpar"=ergrell)
}
