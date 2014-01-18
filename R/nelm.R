nelm <-
  function(reshOBJ, etastart="aut", ctrl=list())
  {
    
    call <-match.call()
    attr(call, "date") <- date()
    
    ## controls
    RA <- qr(reshOBJ$Qmat)$rank
    if(RA < ncol(reshOBJ$Q))
    {
      stop(paste("Q matrix has only rank",RA," - but must have full column rank! \n"))  
    }
    
    #########################################
    ######### USER CONTROLS #################
    #########################################
    
    cont <- list(nodes=14, absrange=5, sigmaest=FALSE, exac=0.00001, EMmax = 500, verbose=TRUE, NRmax=20, NRexac=0.01, dooptim=FALSE, centBETA=FALSE, centALPHA=FALSE,Clist=NA, nonpar=FALSE)
    
    user_ctrlI <- match(names(ctrl),names(cont))
    if(any(is.na(user_ctrlI)))
    {
      notex <- names(ctrl)[is.na(user_ctrlI)]
      warning("The following options in ctrl do not exist: ", paste(notex,collapse=", "))
      ctrl       <- ctrl[!is.na(user_ctrlI)]
      user_ctrlI <- user_ctrlI[!is.na(user_ctrlI)]
    }
    
    cont[user_ctrlI] <- ctrl
    
    #---------------------------------------------------------------------------
    
    
    ## generate starting values
    startOBJ <- startV_nlmMG(reshOBJ=reshOBJ,etastart=etastart,Clist=cont$Clist)
    ##generate quadrature nodes and weights
    quads <- quadIT(nodes=cont$nodes,absrange=cont$absrange,ngr=nlevels(reshOBJ$gr))
    
    
    OLD     <- 0
    ZAEHL   <- 1
    PARS    <- startOBJ$ulstv
    NLev    <- nlevels(reshOBJ$gr)
    ESTlist <- list()
    mueERG <- NA
    
    
    
    
    if(cont$dooptim)
    {
      
      repeat
      {
        if(cont$verbose){cat("\r Estep:",ZAEHL,"| Mstep:", ZAEHL -1,"\r")}
        
        #E ****
        erg_estep <- Enlm(PARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,PREVinp=mueERG)
        
        #M ****
        oerg <- optim(par=PARS,fn=ZFnlm,gr=de1nlm,erg_estep=erg_estep, startOBJ=startOBJ, reshOBJ=reshOBJ, quads=quads, control=list(fnscale=-1,maxit=50),method="BFGS",hessian=TRUE)
        
        
        if(NLev > 1)
        {
          # estimate mu and sigma

          mueERG <- mueNLM(PARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,sigmaest=cont$sigmaest)
          # new quads
          quads <- quadIT(nodes=cont$nodes,absrange=cont$absrange,ngr=NLev,mu=mueERG$mean_est,sigma=mueERG$sig_est)
        }
        
        if(abs(OLD - abs(oerg$value)) <= cont$exac & ZAEHL > 10 | ZAEHL > cont$EMmax)
        {
          mueERG <- mueNLM(PARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,sigmaest=cont$sigmaest,endest=TRUE)
          
          ESTlist[[1]]   <- oerg$par
          ESTlist[[2]]   <- erg_estep
          ESTlist[[3]]   <- oerg
          ESTlist[[4]]   <- ZAEHL
          ESTlist[[5]]   <- mueERG
          ESTlist[[6]]   <- quads
          ESTlist[[7]]   <- startOBJ
          ESTlist[[8]]   <- cont
          
          names(ESTlist) <- c("etapar","last_estep","last_mstep" ,"n_steps","erg_distr","QUAD","starting_values","ctrl")
          break
        }
        
        OLD <- abs(oerg$value)
        PARS <- oerg$par
        ZAEHL <- ZAEHL + 1
      }
      
      
      
      
    } else {
      
      
      
      ## EM Algorithm
      repeat
      {
        if(cont$verbose & (ZAEHL == 1 | NLev <= 1)){cat("\r Estep:",ZAEHL,"| Mstep:", ZAEHL -1,"\r")}
        
        #E ****
        erg_estep <- Enlm(PARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,PREVinp=mueERG,nonpar=cont$nonpar)
        
        if(cont$nonpar) # nonpar estimation in any case (reference group is standardized 0,1)
        {
          quads <- quadIT(nodes=cont$nodes,absrange=cont$absrange,ngr=NLev, ergE=erg_estep)  
        }
        
        #################################
        ## Newton Raphson Procedure  ####
        #################################
        if(cont$verbose){cat("\r Estep:",ZAEHL,"| Mstep:", ZAEHL,"\r")}
        mPARS <- PARS
        
        for(i in 1:cont$NRmax)
        {
          fir1 <- de1nlm(mPARS,erg_estep=erg_estep,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads)
          sec2 <- de2nlm(mPARS,erg_estep=erg_estep,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads)
          
          newP <- mPARS - as.vector(fir1 %*% solve(sec2))
          
          if(sum(abs(newP - mPARS)) < cont$NRexac)
          {
            mPARS <- newP
            break 
          } else {
            mPARS <- newP
          }
        }
        
        

        
        if(sum(abs(mPARS - PARS)) <= cont$exac & ZAEHL > 10 | ZAEHL >= cont$EMmax)
        {
          if(cont$verbose){cat("\r Estep:",ZAEHL,"| Mstep:", ZAEHL," --> estimating EAP & co \r")}
         
          mueERG <- mueNLM(mPARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,sigmaest=cont$sigmaest,endest=TRUE)
          #quads <- quadIT(nodes=cont$nodes,absrange=cont$absrange,ngr=NLev,mu=mueERG$mean_est,sigma=mueERG$sig_est)
          
          ###### LIKELIHOOD BERECHNEN #######################
          value <- ZFnlm(mPARS,erg_estep=erg_estep,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads)
          ###################################################
          
          ESTlist[[1]]   <- mPARS
          ESTlist[[2]]   <- erg_estep
          ESTlist[[3]]   <- list(value=value, hessian=sec2)
          ESTlist[[4]]   <- ZAEHL
          ESTlist[[5]]   <- mueERG
          ESTlist[[6]]   <- quads
          ESTlist[[7]]   <- startOBJ
          ESTlist[[8]]   <- cont
          
          names(ESTlist) <- c("etapar","last_estep","last_mstep" ,"n_steps","erg_distr","QUAD","starting_values", "ctrl")
          break
        }
        
        
        if(NLev > 1 & !cont$nonpar) # more than 1 group and NOT nonpar estimation
        {
          if(cont$verbose & ZAEHL >= 1){cat("\r Estep:",ZAEHL+1,"| Mstep:", ZAEHL,"\r")}  
          # estimate mu and sigma
          #mueERG <- mueNLM(PARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,sigmaest=cont$sigmaest)
          mueERG <- mueNLM(mPARS,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads,sigmaest=cont$sigmaest)
          # new quads
          quads <- quadIT(nodes=cont$nodes,absrange=cont$absrange,ngr=NLev,mu=mueERG$mean_est,sigma=mueERG$sig_est)
        }
        
        ###!
        value <- ZFnlm(mPARS,erg_estep=erg_estep,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads)
        
        cat("logL: ",value,"\n")
        
        PARS <- mPARS
        ZAEHL <- ZAEHL + 1
        
        
      }
      
    }

    ## Person Parameters
    ESTlist$EAPs <- mueERG$thetas
    # center the parameters
    parcent <- Cnlm(reshOBJ=reshOBJ,ESTlist=ESTlist, centBETA=cont$centBETA, centALPHA=cont$centALPHA, startOBJ = startOBJ)
    
    
    ESTlist$ZLpar <- parcent
    
    #SE estimation ----------------
    
    if(all(!is.na(startOBJ$setC)))
    {
      
      comphess <- diag(reshOBJ$Qmat[,-startOBJ$setC$whichetas] %*% solve(ESTlist$last_mstep$hessian) %*% t(reshOBJ$Qmat[,-startOBJ$setC$whichetas]))
      comphesq <- sqrt(comphess*(-1))
      
      ## hier noch die constanten auf NA setzen!

      comphesq[rowSums(cbind(reshOBJ$Qmat[,startOBJ$setC$whichetas]))  > 0] <- NA
      
      
    } else 
    {
      comphess <- diag(reshOBJ$Qmat %*% solve(ESTlist$last_mstep$hessian) %*% t(reshOBJ$Qmat))
      comphesq <- sqrt(comphess*(-1))
      #
    }
    
    
    #comphess <- diag(reshOBJ$Qmat %*% solve(ESTlist$last_mstep$hessian) %*% t(reshOBJ$Qmat))
    #comphesq <- sqrt(comphess*(-1))
    #
    notest <- which(rowSums(reshOBJ$Qmat) == 0)
    comphesq[notest] <- NA
    attr(comphesq,"listform")  <- relist(comphesq,startOBJ$stwm1)
    ESTlist$SE <- comphesq
#)
    # also save the reshOBJ - and of course the queen
    ESTlist$reshOBJ   <- reshOBJ

    # -----------------------------
    # adding category information
    # -----------------------------
    ESTlist$Catinf <- Infnelm(ESTlist)
    
    ESTlist$call <- call
    
    class(ESTlist) <- "nelm"
    ## write out
    return(ESTlist)  
  }
