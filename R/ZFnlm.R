ZFnlm <-
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
  
  
  opp   <- as.vector(Q %*% Ulstv)
  
  relstv <- relist(opp,SKEL)
  
 
  ALL_G <- mapply(function(stvl,ql,RI)
        {
          
         riq_quer_mat <- sapply(RI$riq_querG,function(x)x)
         f_iq         <- sapply(RI$riqv_querG,colSums) + riq_quer_mat # nodes x items
         
         ZQsternl <- lapply(stvl,function(pitem)
         {   
           # 2PL Part
           # ----------------------
           Km      <- matrix(c(rep(1,length(ql$nodes)),ql$nodes),ncol=2)
           
           woba <- grep("(beta|alpha)",names(pitem))
           abpar <- pitem[woba]
           
           solit <- twoplpart(Km=Km, abpar=abpar)

           dosolit <- 1-solit
           
           # NRM Part
           # ----------------------
           pitemNRM <- pitem[-woba]
           
           ZQstern <- coP_nrm(pitemNRM,Km)
           
           ZQstern_all <- cbind(solit,ZQstern)
           ZQstern_all
         }) 
         
         
         logPiXq <- sapply(ZQsternl,function(x)log(x[,1]))
         T1      <- sum(logPiXq  * riq_quer_mat)
         
         sumriqv_quer    <- sapply(RI$riqv_querG,function(x) colSums(x))
         log_EMINUS_PiXq <- sapply(ZQsternl,function(x)log(1 - x[,1]))
         T2              <- sum(log_EMINUS_PiXq * sumriqv_quer)
         
         logPivuis0 <- lapply(ZQsternl,function(x)log(x[,-1])) 
         T3         <- sum(mapply(function(A,B){sum(t(A) * B)},A=RI$riqv_querG,B=logPivuis0))
         
         allin <- T1 + T2 + T3
         allin

          },stvl=relstv,ql=quads,RI=erg_estep,SIMPLIFY = FALSE)
  
 endsum <- sum(sapply(ALL_G,sum))
  return(endsum)
 
}
