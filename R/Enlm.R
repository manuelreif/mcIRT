Enlm <-
function(Ulstv,reshOBJ,startOBJ,quads,PREVinp,nonpar)
{ 
  # new
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
    
  Km      <- matrix(c(rep(1,length(ql$nodes)),ql$nodes),ncol=2) 
    
  nrmez <- mapply(function(pitem,daten)
      {
        # 2PL Part
        # ----------------------
        
        woba <- grep("(beta|alpha)",names(pitem))
        abpar <- pitem[woba]
        
        solit <- twoplpart(Km=Km, abpar=abpar)
    
        
        dosolit <- 1-solit
        tplpart <- cbind(solit,dosolit)
        
        # NRM Part
        # ----------------------
        pitemNRM <- pitem[-woba]
        
        ZQstern <- coP_nrm(pitemNRM,Km)
        
        ZQstern_nlm <- ZQstern * as.vector(dosolit)
        
        # now the data
        # ---------------------
        ZQstern_all <- cbind(solit,ZQstern_nlm)
        
        ZQdstern <- ZQstern_all %*% t(daten)

        return(ZQdstern=ZQdstern)
      },pitem=stvl,daten=d1,SIMPLIFY = FALSE)
  
  },stvl=relstv, ql=quads, levs=levels(reshOBJ$gr), d1=dAtA, SIMPLIFY = FALSE)
  
  
  # ----------------------------------------------------------- 
  
  riq_querA <- mapply(function(levs, zqgroup, ql, d1uc)
  {
    nrme2   <- array(unlist(zqgroup),c(dim(zqgroup[[1]]),length(zqgroup)))
    
    bullet_with_butterfly_wings <- apply(nrme2,1:2,function(x) prod(x,na.rm=T))
    # nodes x persons
    
    #anzahlnodes <- 1:dim(nrme2)[1]
    LjmalA <- bullet_with_butterfly_wings * ql$weights
    Pji_schlange  <- colSums(LjmalA) 
    
    ################################################
    #-------------------------- riqv_quer -----------
    ################################################
    riqv_1teil  <- t(LjmalA) / Pji_schlange
    
    
    riq_querG <- lapply(d1uc,function(x)
    {
      t(x[,1]) %*% riqv_1teil
    })
    
    riqv_querG <- lapply(d1uc,function(x)
    {
      t(x[,-1]) %*% riqv_1teil
    })
    
    #list(riq_querG=riq_querG, riqv_querG=riqv_querG)
    
    if(nonpar)
    {
      return(list(riq_querG=riq_querG, riqv_querG=riqv_querG,fquer=riqv_1teil)) # fquer added for nonpar distr estimation
    } else {
            return(list(riq_querG=riq_querG, riqv_querG=riqv_querG))  
           }
    
    
  },levs=levels(reshOBJ$gr), zqgroup=nrme1, ql=quads, d1uc=datuc, SIMPLIFY = FALSE)
  ############xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  }  else 
  {
    
    riq_querA <- mapply(function(levnr, ql, d1uc, MER)
    {
      
      ################################################
      #-------------------------- riqv_quer -----------
      ################################################
      riqv_1teil  <- t(MER$LjmalA) / MER$Pji_schlange   
      
      
      riq_querG <- lapply(d1uc,function(x)
      {
        t(x[,1]) %*% riqv_1teil
      })
      
      riqv_querG <- lapply(d1uc,function(x)
      {
        t(x[,-1]) %*% riqv_1teil
      })
      
      
      if(nonpar)
      {
       return(list(riq_querG=riq_querG, riqv_querG=riqv_querG,fquer=riqv_1teil)) # fquer added for nonpar distr estimation
      } else {
             return(list(riq_querG=riq_querG, riqv_querG=riqv_querG))  
             }
      
      #return(list(riq_querG=riq_querG,riqv_querG=riqv_querG,fquer=riqv_1teil))
    },levnr=1:length(levels(reshOBJ$gr)),ql=quads, d1uc=datuc,MER=PREVinp$mue_hat_g,SIMPLIFY = FALSE)
    
    
  }
  
  return(riq_querA=riq_querA)

}
