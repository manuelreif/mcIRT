Infnrm <- function(ESTlist, fromto=c(-5,5), gran=200)
{
  # ESTlist --> internal ESTlist created by nrm. must contain the estimated values which stem from the EM procedure as well as the centered category parameters
  
  #SKEL  <- ESTlist$starting_values$stwm1
  #Q     <- ESTlist$reshOBJ$Qmat
  #opp   <- as.vector(Q %*%  ESTlist$etapar)
  
  # wahrscheinlich ist es besser die ZLpar list zu nehmen.
  
  #relstv <- relist(opp,SKEL)
  relstv  <- ESTlist$ZLpar
  thetas <- seq(fromto[1], fromto[2], length.out=gran)
  
#   levs=levels(ESTlist$reshOBJ$gr)[1]
#   stvl=relstv[[1]]
#   ql=ESTlist$QUAD[[1]]
  
  catinfG <- mapply(function(levs,stvl)
  { # loops all groups
    
    # pitem = stvl[[1]]
    catinfI <- mapply(function(pitem)
    { # loops all items
      
      Km  <- matrix(c(rep(1,length(thetas)), thetas),ncol=2)
      LAM <- matrix(pitem,nrow=2,byrow=T)
      
      Z <- Km %*% LAM
      ez <- exp(Z)
      ezrs <- rowSums(ez)        
      ZQstern <- ez / ezrs
      
      # ZQstern
      # PROB = ZQstern
      # whichI = itemnr
      # TT = gibts nicht - es kann erst ausserhalb der schleife mit der Q Matrix multipliziert werden!
      
      #TT <- PROB[[2]]
      #     
      # zei <- ZQstern[1,]
      # zei <- 1
      LAMs <- pitem[(length(pitem)/2 + 1):length(pitem)]
      
      W_g <- sapply(1:nrow(ZQstern),function(zei) # geht die nodes durch
      {
        Zqrow <- ZQstern[zei,]
        z     <- thetas[zei]
        Pqrep <- matrix(-Zqrow, length(Zqrow), length(Zqrow)) 
        diag(Pqrep) <- 1-Zqrow 
        Pdi <- diag(Zqrow)
        Wi  <- Pqrep %*% Pdi
        as.vector(LAMs %*% Wi %*% LAMs) * ZQstern[zei,]
        
      }) 
      
     t(W_g) # Category Information
      
    },pitem=stvl,  SIMPLIFY = F) ### was T before
    
    catinfI
  },levs=levels(ESTlist$reshOBJ$gr), stvl=relstv ,SIMPLIFY = FALSE)
  
  
  # category informations - for different thetas for each group
  class(catinfG) <- "infnrm"
  
  
  return(list(catinfG=catinfG, thetas=thetas))
  
}  



