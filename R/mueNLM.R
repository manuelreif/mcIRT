mueNLM <-
function(Ulstv,reshOBJ,startOBJ,quads,sigmaest=FALSE,endest=FALSE)
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
  
  
  
  nrme1 <- mapply(function(stvl,ql,levs, d1)
  { # loops all groups
    
  Km      <- matrix(c(rep(1,length(ql$nodes)),ql$nodes),ncol=2)
    
    nrmez <- mapply(function(pitem,daten) # rattert pro item durch
    {
      # 2PL Part
      # ----------------------
      
      woba  <- grep("(beta|alpha)",names(pitem))
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
      # missing values remain as NA in the Data-Structure
      
      return(ZQdstern=ZQdstern)
    },pitem=stvl,daten=d1,SIMPLIFY = FALSE)
    
  },stvl=relstv, ql=quads, levs=levels(reshOBJ$gr), d1=dAtA,SIMPLIFY = FALSE)
  
  
  
## ------------  

  
  
  mue_hat_g <- mapply(function(levs, zqgroup, ql)
  {

    nrme2   <- array(unlist(zqgroup),c(dim(zqgroup[[1]]),length(zqgroup)))
    
    bullet_with_butterfly_wings <- apply(nrme2,1:2,function(x) prod(x,na.rm=T))
    # nodes x persons
    

    LjmalA <- bullet_with_butterfly_wings * ql$weights
    Pji_schlange  <- colSums(LjmalA) 
    

    
    LjmalAmalN <- bullet_with_butterfly_wings * ql$weights * ql$nodes
    Xquer_insum <- colSums(LjmalAmalN)
    
    Xqi_quer  <- Xquer_insum / Pji_schlange   
    Xqi_qq    <- mean(Xqi_quer)
    
    # SIGMA
    if(sigmaest & length(levels(reshOBJ$gr)) > 1)
    {
      
      si_term1 <- t(sapply(ql$nodes,function(NOO) (NOO - Xqi_quer)^2)) 
      
      sigma_xugiZ <- colSums(si_term1 * bullet_with_butterfly_wings * ql$weights)
      sigma_xugi  <- sigma_xugiZ / Pji_schlange
      
      
      Xq_min_muHat <- (Xqi_quer - Xqi_qq)^2
      
      SIGhatsq     <- mean(sigma_xugi + Xq_min_muHat)#
      
      mesi <- c(Xqi_qq, sqrt(SIGhatsq))
      
    } else {
      mesi <- c(Xqi_qq, 1)  
    }
    
    
    ####################### SE ##################

    if(endest & sigmaest) 
    {
      LdurchsumL <- t(LjmalA) / Pji_schlange 
      sumOnodes <- colSums(t(LdurchsumL) * (ql$nodes - Xqi_qq)/mesi[2])^2
      muemue    <- sum(sumOnodes)
      
      sumOsignodes <- colSums(t(LdurchsumL) * ((ql$nodes - Xqi_qq)^2 - mesi[2]) / (2*mesi[2]^2))^2
      sigsig <- sum(sumOsignodes)
      
      muesig <- colSums(t(LdurchsumL) * (ql$nodes - Xqi_qq)/mesi[2]) %*% colSums(t(LdurchsumL) * ((ql$nodes - Xqi_qq)^2 - mesi[2]) / (2*mesi[2]^2))
      
      Infmat <- matrix(c(muemue,muesig,muesig,sigsig),2,2)
      serrors <- sqrt(diag(solve(Infmat)))
      names(serrors) <- c("SE_mue","SE_sigma")
    } else if(endest){ # if EM is converged and no sigma estimation took place
      LdurchsumL <- t(LjmalA) / Pji_schlange 
      sumOnodes <- colSums(t(LdurchsumL) * (ql$nodes - Xqi_qq)/mesi[2])^2
      muemue    <- sum(sumOnodes)
      serrors     <- c(sqrt(1/muemue),NA)
      names(serrors) <- c("SE_mue","SE_sigma")
    } else  {
      serrors <- NA
    }
    
    
    names(mesi) <- paste(c("mean","var"),"_",levs,sep="")
    
    #
    list(mesi=mesi,LjmalA=LjmalA,Pji_schlange=Pji_schlange,serrors=serrors)
    

  },levs=levels(reshOBJ$gr), zqgroup=nrme1, ql=quads, SIMPLIFY = FALSE)
  ############xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  # centering 
  getmean <- sapply(mue_hat_g,function(x)x$mesi[1])
  getvar <- sapply(mue_hat_g,function(x)x$mesi[2])
  ##
  mean_est <- (getmean - getmean[1])
  sig_est  <- getvar + (1-getvar[1])
  
  #######
  if(endest){
    errmat <- sapply(mue_hat_g,function(x) x$serrors)
  } else {errmat <- NA}
  
  return(list(mean_est=mean_est , sig_est=sig_est, mue_hat_g=mue_hat_g, errmat=errmat))
}
