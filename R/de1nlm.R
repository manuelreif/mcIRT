de1nlm <-
function(Ulstv,erg_estep,startOBJ,reshOBJ,quads)
{
  
  SKEL  <- startOBJ$stwm1
  Q     <- reshOBJ$Qmat
  
  
  if(all(!is.na(startOBJ$setC)))
  {
    
    bigv <- vector(mode="numeric",length=ncol(Q))
    
    bigv[-startOBJ$setC$whichetas] <- Ulstv
    bigv[startOBJ$setC$whichetas]  <- startOBJ$setC$whichconstant
    
    Ulstv <- bigv
  }
  
  
  opp    <- as.vector(Q %*% Ulstv)
  relstv <- relist(opp,SKEL)
  
  
  ALL1st <- mapply(function(stvl,ql,RI)
  { # loops all groups
    
  Km      <- matrix(c(rep(1,length(ql$nodes)),ql$nodes),ncol=2) #!
    
   
      deriv4IT <- mapply(function(pitem,II) # rattert pro item durch
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
          ZQstern_all <- cbind(solit,ZQstern)
          
    
          
          # 2pl part
          # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          fique0 <- sapply(RI$riqv_querG,colSums)
          riq_quer_mat <- sapply(RI$riq_querG,function(x)x)
          f_iq         <- fique0 + riq_quer_mat # nodes x items
          
            # II denotes the number of the current item
            commonterm <- as.vector(RI$riq_querG[[II]] - ZQstern_all[,1] * f_iq[,II]) ## II brauchst du schon
            
            betaP <- sum(commonterm)
            names(betaP) <- paste("beta_IT",II,sep="")
            
            alphaP <- sum(commonterm * ql$nodes)
            names(alphaP) <- paste("alpha_IT",II,sep="")
          
          # NRM part
          # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
          #fique0 <- sapply(RI$riqv_querG,colSums)
          

            fqomega <- fique0[,II] * ZQstern_all[,-1]
            
            Rfo     <- t(RI$riqv_querG[[II]]) - fqomega
          
            GAMderiv <- colSums(Rfo)        
            XIderiv  <- colSums(Rfo *ql$nodes)

            
            c(betaP,GAMderiv,alphaP,XIderiv)
            

        },pitem=stvl,II=1:length(RI$riq_querG),SIMPLIFY = FALSE)
      
      deriv4IT
        
  },stvl=relstv,ql=quads,RI=erg_estep,SIMPLIFY = FALSE)


  deriv <- unlist(ALL1st)
  
  derivV <- as.vector(deriv %*% Q)
  names(derivV) <- colnames(Q)
  
  if(all(!is.na(startOBJ$setC)))
  {
    derivV <- derivV[-startOBJ$setC$whichetas]
  }
  
  
  
  derivV

}
