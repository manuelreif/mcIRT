Enrm <-
function(Ulstv,reshOBJ,startOBJ,quads,PREVinp)
{

  
  dAtA    <- reshOBJ$recm
  datuc   <- reshOBJ$d1uc
  
  SKEL  <- startOBJ$stwm1
  Q     <- reshOBJ$Qmat
  
  # in case some paraters are set to a certain value
  if(all(!is.na(startOBJ$setC)))
    {
      
    bigv <- vector(mode="numeric",length=ncol(Q))
      
    bigv[-startOBJ$setC$whichetas] <- Ulstv
    bigv[startOBJ$setC$whichetas]  <- startOBJ$setC$whichconstant
    
    Ulstv <- bigv
    }
  
  
  opp    <- as.vector(Q %*% Ulstv)
  relstv <- relist(opp,SKEL)
  
if(all(is.na(PREVinp)))
{ 
  
  nrme1 <- mapply(function(stvl,ql,levs, d1)
  { # loops all groups
    
    Km  <- matrix(c(rep(1,length(ql$nodes)),ql$nodes),ncol=2)

    nrmez <- mapply(function(pitem,daten) 
    { # loops all items
      

      LAM <- matrix(pitem,nrow=2,byrow=T)
      
      ez      <- exp(Km %*% LAM)
      ezrs    <- rowSums(ez)        
      ZQstern <- ez / ezrs
      
      ZQdstern <- ZQstern %*% t(daten)
      
      return(ZQdstern)
    },pitem=stvl, daten=d1, SIMPLIFY = FALSE)


  },stvl=relstv, ql=quads, levs=levels(reshOBJ$gr), d1=dAtA, SIMPLIFY = FALSE)
  
  
  riqv_quer <- mapply(function(levs, zqgroup, ql, d1uc)
  {
    nrme2   <- array(unlist(zqgroup),c(dim(zqgroup[[1]]),length(zqgroup))) #recreate it as 3d array
    # nodes x persons(groupX) x items
    
    bullet_with_butterfly_wings <- apply(nrme2,1:2,function(x) prod(x,na.rm=TRUE))
    # nodes x persons
    
    anzahlnodes <- 1:dim(nrme2)[1]
    LjmalA <- bullet_with_butterfly_wings * ql$weights
    Pji_schlange  <- colSums(LjmalA) 
    
    ################################################
    #-------------------------- riqv_quer -----------
    ################################################
    riqv_1teil  <- t(LjmalA) / Pji_schlange    
    riqv_querG <- lapply(d1uc,function(x)
    {
      t(x) %*% riqv_1teil
    })
    
    riqv_querG
  },levs=levels(reshOBJ$gr), zqgroup=nrme1, ql=quads, d1uc=datuc, SIMPLIFY = FALSE)
  
} else { 
 
        riqv_quer <- mapply(function(levs, d1uc, ql, MER)
        {
          ################################################
          #-------------------------- riqv_quer -----------
          ################################################
          
          riqv_1teil  <- t(MER$LjmalA) / MER$Pji_schlange   
          riqv_querG <- lapply(d1uc,function(x)
          {
            t(x) %*% riqv_1teil
          })
          
          riqv_querG
        },levs=levels(reshOBJ$gr), d1uc=datuc, ql=quads, MER=PREVinp$mue_hat_g, SIMPLIFY = FALSE)  
        

      }
  
  return(riqv_quer=riqv_quer)
}
