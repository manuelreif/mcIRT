Infnelm <- function(ESTlist, fromto=c(-5,5), gran=200)
{
  # ESTlist --> internal ESTlist created by nrm. must contain the estimated values which stem from the EM procedure as well as the centered category parameters
  
  #SKEL  <- ESTlist$starting_values$stwm1
  #Q     <- ESTlist$reshOBJ$Qmat
  #opp   <- as.vector(Q %*%  ESTlist$etapar)
  
  # wahrscheinlich ist es besser die ZLpar list zu nehmen.
  
  #relstv <- relist(opp,SKEL)
  
  #ANSTATT "relstv"
  ALPHAS  <- ESTlist$ZLpar$alpha   # alphas für alle items
  BETAS   <- ESTlist$ZLpar$beta    # betas für alle items
  NRMS    <- ESTlist$ZLpar$nrmpar  # nrm parameters für alle items
  
  # thetas - x-axis for the cat information plot
  thetas <- seq(fromto[1], fromto[2], length.out=gran)
  
  #   levs=levels(ESTlist$reshOBJ$gr)[1]
#   alphaG=ALPHAS[[1]]
#   betaG = BETAS[[1]]
#   nrmG = NRMS[[1]] 
  
  
  catinfG <- mapply(function(levs, alphaG, betaG, nrmG)
  { # loops all groups
    
    
#     alI  = alphaG[1]
#     betI = betaG[1]
#     nrmI = nrmG[[1]]
    
    catinfI <- mapply(function(alI, betI, nrmI)
    { # loops all items
      
      pitem <- c(nrmI$zetas,nrmI$lambdas)
      Km  <- matrix(c(rep(1,length(thetas)), thetas),ncol=2)
      LAM <- matrix(pitem,nrow=2,byrow=T)
      
      Z <- Km %*% LAM
      ez <- exp(Z)
      ezrs <- rowSums(ez)        
      ZQstern <- ez / ezrs
      
      # 2PL
      abpar <- c(betI,alI)
      P2pl <- twoplpart(Km=Km, abpar=abpar)
      Q2pl <- (1-P2pl)
      
      TwoPLInf <- alI^2 * P2pl * P2pl * Q2pl # sihe Psychometrika Artikel S.462
      #####
      #### Normally the item information function is:
      # alpha^2 * P * Q
      # but for the Nested logit model it is alpha^2 * P * P * Q --> because you are modeling the distractors as well.
      # A similar thing is true for the nrm - part. as we know, the nrm is modeled in case the correct answer was not found --> so it depends on 1-P !
      # --> this can be seen in line labels with: # ***
      
      
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
        as.vector(LAMs %*% Wi %*% LAMs) * Q2pl[zei] * ZQstern[zei,] # ***
        
      }) 
      

      
      
      cbind(TwoPLInf,t(W_g)) # Category Informations incl 2PL part
      
    },alI = alphaG, betI = betaG, nrmI = nrmG,  SIMPLIFY = F) ### was T before
    
    
    catinfI
    
  },levs=levels(ESTlist$reshOBJ$gr), alphaG=ALPHAS ,betaG = BETAS, nrmG = NRMS ,SIMPLIFY = FALSE)
  
  
  TIFall <- lapply(catinfG,function(GRs)
              {
              apply(simplify2array(GRs, higher=TRUE),1,sum)
              })

  
  # category informations - for different thetas for each group
  class(catinfG) <- "infnlm"
  
  
  return(list(catinfG=catinfG, thetas=thetas, TestInfGROUPS=TIFall))
  
}  


