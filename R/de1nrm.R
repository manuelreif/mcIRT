de1nrm <-
function(Ulstv,riqv_quer,startOBJ,reshOBJ,quads)
{
  
  SKEL  <- startOBJ$stwm1
  Q     <- reshOBJ$Qmat
  opp   <- as.vector(Q %*% Ulstv)
  
  relstv <- relist(opp,SKEL)
  
  ## 1
  fiq <- lapply(riqv_quer,function(X) sapply(X,function(newx)colSums(newx)))
  
#   
  occ <- mapply(function(stvl,ql,levs,RI,FI)
  { # loops all groups
    
    Km  <- matrix(c(rep(1,length(ql$nodes)),ql$nodes),ncol=2)
    
    nrmez <- mapply(function(pitem,rii,itnr)
    { # loops all items
      
      LAM <- matrix(pitem,nrow=2,byrow=T)
      
      Z <- Km %*% LAM
      ez <- exp(Z)
      ezrs <- rowSums(ez)        
      ZQstern <- ez / ezrs
      
      fqomega <- FI[,itnr] * ZQstern
      Rfo     <- t(rii) - fqomega
      
      GAMderiv <- colSums(Rfo)        
      XIderiv  <- colSums(Rfo *ql$nodes)
      
      c(GAMderiv,XIderiv)
      
    },pitem=stvl,rii=RI,itnr=1:ncol(FI),SIMPLIFY = F) 
    
    nrmez
    
  },levs=levels(reshOBJ$gr),stvl=relstv,ql=quads,RI=riqv_quer,FI=fiq,SIMPLIFY = FALSE)
  
  
  deriv <- unlist(occ)
  
  derivV <- as.vector(deriv %*% Q)
  names(derivV) <- colnames(Q)
  derivV
}
