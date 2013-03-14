grDM <-
function(aDD,gr,design,TYPE=TYPE)
  
# design explained:
# the input for the design argument could be a character like "nodif" - if so, the Q matrix is built under the assumption that there is no dif to be measured between groups
  # on the other hand the input for design could also be a list.
   # the first list element must be: $zeta ; the second must be $lambda and refers to the two parameter types in the nominal response model
    # in each list slot - a matrix of  NUMBEROFGROUPS x NUMBEROFITEMS is expected which contains information about which item has to be estimated in which group serperately.
      # the first row of the matrix (in each of the 2 slots) has to be 1 only - which means, that in the first group every parameter is estimated
        # the 2nd row determines whether the parameters (zeta/lambda) are estimated new for this group, or are the same as in the first group. when they are supposed to be the same, it is reasonable to estimate them together with the first group --> fill in 1. if you like to estimate the parameter seperately for this group, fill in the number of the grou --> 2. so for 3 groups the matrix (eg for the zetas) could look like this:

# bspmatrix <- matrix(c(1,1,1,1,1,1,1,1,2,2,1,1,1,2,3),nrow=3,byrow=T)
# 
# [,1] [,2] [,3] [,4] [,5]
# [1,]    1    1    1    1    1
# [2,]    1    1    1    2    2
# [3,]    1    1    1    2    3

# so this matrix means, that the zetas in item 1 - 3 are estimated together in all groups,
# zetas in item 4 are estimated in group 1 and seperately in group 2 and 3
# and the zetas in item 5 are estimated seperately in each group - so there will be 3 different estimates of the zetwas for this item !
  
{
  
  #ctrl
  
  if(any(is.na(gr))){stop("Missing values inside the grouping vector are not allowed")}
  
  if(TYPE=="NRM")
  {
  
  if(nlevels(gr) == 1) # if no groups
  {
    cats <- sapply(aDD,function(x)x$anz_cat)
    nro  <- sum(cats*2)  
    Qo    <- matrix(0,nro,nro)
    diag(Qo) <- 1
    
    leaveOUT <- cumsum(rep(cats,each=2))
    Q <- Qo[,-leaveOUT]
    
    #naming
    prae <- rep(paste("I",1:length(cats),sep=""),cats*2)
    app1 <- unlist(sapply(cats,function(AA) paste(rep(c("zeta","lam"),each=AA),rep(1:AA,2),sep="")))
    rwn  <- paste(prae,app1,sep="")
    rownames(Q) <- rwn
    colnames(Q) <- paste("eta",1:ncol(Q),sep="")
    # DONE
    
  } else if(nlevels(gr) > 1) # if there are any groups (> than 1 group)!
  {
    numbcat <- length(levels(gr))
    cats    <- sapply(aDD,function(x)x$anz_cat)
    nro     <- sum(cats*2)  
    Qbbo     <- matrix(0,nro,nro) 
    diag(Qbbo) <- 1
    mult1   <- matrix(0,numbcat,numbcat)
    
    leaveOUT <- cumsum(rep(cats,each=2))
    Qbb <- Qbbo[,-leaveOUT]
    
    # names with multiple groups
    prae1 <- rep(rep(paste("I",1:length(cats),sep=""),cats*2),numbcat)
    prae <- paste(rep(paste("G",1:numbcat,sep=""),each=nro),prae1,sep="")
    app1 <- rep(unlist(sapply(cats,function(AA) paste(rep(c("zeta","lam"),each=AA),rep(1:AA,2),sep=""))),numbcat)
    rwn  <- paste(prae,app1,sep="")
    
    ##################
    #  DIF ? #########
    ##################
    
    if(all(design == "nodif"))
    {
      mult1[,1] <- 1
      Q <- mult1 %x% Qbb
      onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
      Q <- Q[,-onlyZ]
      rownames(Q) <- rwn
      colnames(Q) <- paste("eta",1:ncol(Q),sep="")
      
      
    } else if(all(design == "dif1"))
    {
      # design with big identity matrix
      diag(mult1) <- 1
      Q <- mult1 %x% Qbb
      
      
      rownames(Q) <- rwn
      colnames(Q) <- paste("eta",1:ncol(Q),sep="")
      
      
    } else if(all(design == "dif2"))
    {
      # design where zetas are estimated within each group - lambdas remaining the same for each group
      mult1[lower.tri(mult1,diag=TRUE)] <- 1
      Q <- mult1 %x% Qbb
      
      whzeta <- grep("zeta",rwn)
      whzeta1 <- whzeta[whzeta > nrow(Qbb)]
      
      Q[whzeta1,1:ncol(Qbb)] <- 0 
      
      whlam <- grep("lam",rwn)
      whlam1 <- whlam[whlam > nro]
      Q[whlam1,-(1:ncol(Qbb))] <- 0
      
      onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
      Q <- Q[,-onlyZ]
      
      rownames(Q) <- rwn
      colnames(Q) <- paste("eta",1:ncol(Q),sep="")
      
      
    } else if(all(design == "dif3"))
    {
      # design where lambda is estimated for each group while the zetas remain the same for all different groups
      mult1[lower.tri(mult1,diag=TRUE)] <- 1
      Q <- mult1 %x% Qbb
      
      whzeta <- grep("lam",rwn)
      whzeta1 <- whzeta[whzeta > nrow(Qbb)]
      
      Q[whzeta1,1:ncol(Qbb)] <- 0 
      
      whlam <- grep("zeta",rwn)
      whlam1 <- whlam[whlam > nro]
      Q[whlam1,-(1:ncol(Qbb))] <- 0
      
      onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
      Q <- Q[,-onlyZ]
      
      rownames(Q) <- rwn
      colnames(Q) <- paste("eta",1:ncol(Q),sep="")
      
    } else if(is.list(design))
    {
      
      # HIER FEHLT NOCH EINE KONTROLLFUNKTION FÜR DIE DESIGNLISTE!!!
      #########
      ##################################
      
      # hier jetzt falls ein eigenes design angegeben wird.
      #mult1[,1] <- 1
      mult1[lower.tri(mult1,diag=TRUE)] <- 1
      
      Q <- mult1 %x% Qbb
      rownames(Q) <- rwn

      
      spaltBEG <- c(0,rep(ncol(Qbb),nlevels(gr)) * 1:nlevels(gr))[-(nlevels(gr)+1)] +1 
      spaltEND <- c(rep(ncol(Qbb),nlevels(gr)) * 1:nlevels(gr))
      #aDD
      zeilBEG <- c(0,rep(nrow(Qbb),nlevels(gr)) * 1:nlevels(gr))[-(nlevels(gr)+1)] +1
      zeilEND <- c(rep(nrow(Qbb),nlevels(gr)) * 1:nlevels(gr))
      
      prmA <- c("zeta","lam")
      

                   # gehe die zeilen durch
                    #begz=zeilBEG[2]
                    #endz=zeilEND[2]
                    #
                    zerg <- mapply(function(begz,endz,galaZ) # speichert alle zeilen
                            {
                            #begsp=spaltBEG[1]
                            #endsp=spaltEND[1]
                            # geht die spalten durch:
                            prot <- mapply(function(begsp,endsp,grunr) # speichert je eine vollständige zeile
                                {
                                
                                Qtemp <- Q[begz:endz,begsp:endsp] #auswahl des quadranten
                                
                              
                              for(EACH in 1:length(design))
                                {
                                gala  <- design[[EACH]]
                                loe   <- which(gala[galaZ,] != grunr)
                                prm   <- prmA[EACH]
                                
                                  if(length(loe) == 0)
                                    {
                                      #Qtemp
                                      next
                                      
                                    } else {
  
                                        iaus <- paste("^.+(",paste("I",loe,collapse="|",sep=""),")",prm,sep="")
                                        wo <- grep(iaus,rownames(Qtemp),value=F,perl=TRUE)
                                        Qtemp[wo,] <- 0

                                            }
                                }
                                Qtemp
                                },begsp=spaltBEG, endsp=spaltEND, grunr=1:nlevels(gr),SIMPLIFY=FALSE)
                      
                      

                            },begz=zeilBEG, endz=zeilEND,galaZ=1:nlevels(gr),SIMPLIFY=FALSE)

      
     zwQ <- lapply(zerg,function(XX)do.call(cbind,XX))
      Q  <- do.call(rbind,zwQ)
      onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
      Q <- Q[,-onlyZ]
      
      rownames(Q) <- rwn
      colnames(Q) <- paste("eta",1:ncol(Q),sep="")
      
    } else {stop("Check your input for argument: 'design' ")}
    
  } 
  
  }
  #Q 
  
  #### necessary changes:
  # 1) at least the rownames of matrix Q have to be changed in the NLM model
  # 2) change DIF settings: DIF means different a,b parameters in groups
  #### DDF means, different zeta, lambdas in groups
  # 3) So the design has to change - there must be another 2 list elements
  #### one of "a" and one for "b"
  ##############################################################################
  
  
  
  # -----------------------------------------------#
  ############## ---->  NLM <---- ##################
  ##################################################
  
  
  # übergangs funktion zum testen:
  # funktioniert jetzt natürlich nur bei einer gruppe (also v.a. bei nodif!)
  
  if(TYPE=="NLM")
  {
    
    if(nlevels(gr) == 1) # if no groups
    {
      cats <- sapply(aDD,function(x)x$anz_cat)
      nro  <- sum(cats*2)  
      Qo    <- matrix(0,nro,nro)
      diag(Qo) <- 1
      
      leaveOUT <- cumsum(rep(cats,each=2))
      Q <- Qo[,-leaveOUT]
      
      #naming
      prae <- rep(paste("I",1:length(cats),sep=""),cats*2)
      app1 <- unlist(sapply(cats,function(AA) paste(rep(c("zeta","lam"),each=AA),rep(0:(AA-1),2),sep=""),simplify=FALSE ))
      app1[grep("zeta0",app1)] <- "beta"
      app1[grep("lam0",app1)]  <- "alpha"
      
      rwn  <- paste(prae,app1,sep="")
      rownames(Q) <- rwn
      colnames(Q) <- paste("eta",1:ncol(Q),sep="")
      # DONE
      
    } else if(nlevels(gr) > 1)
    {
      
      numbcat <- nlevels(gr)
      cats    <- sapply(aDD,function(x)x$anz_cat)
      nro     <- sum(cats*2)  
      Qbbo     <- matrix(0,nro,nro) 
      diag(Qbbo) <- 1
      mult1   <- matrix(0,numbcat,numbcat)
      
      leaveOUT <- cumsum(rep(cats,each=2))
      Qbb <- Qbbo[,-leaveOUT]
      
      # names with multiple groups
      prae1 <- rep(rep(paste("I",1:length(cats),sep=""),cats*2),numbcat)
      prae <- paste(rep(paste("G",1:numbcat,sep=""),each=nro),prae1,sep="")
      
      app1 <- rep(unlist(sapply(cats,function(AA) paste(rep(c("zeta","lam"),each=AA),rep(0:(AA-1),2),sep=""))),numbcat)
      app1[grep("zeta0",app1)] <- "beta"
      app1[grep("lam0",app1)]  <- "alpha"
      
      rwn  <- paste(prae,app1,sep="")


      ##################
      #  DIF ? #########
      ##################   
      
      if(all(design == "nodif"))
      {
        mult1[,1] <- 1
        Q <- mult1 %x% Qbb
        onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
        Q <- Q[,-onlyZ]
        rownames(Q) <- rwn
        colnames(Q) <- paste("eta",1:ncol(Q),sep="")
        
        
      } else if(is.list(design))
          {
          
          # hier jetzt falls ein eigenes design angegeben wird.
          #mult1[,1] <- 1
          mult1[lower.tri(mult1,diag=TRUE)] <- 1
          
          Q <- mult1 %x% Qbb
          rownames(Q) <- rwn
          
          
          spaltBEG <- c(0,rep(ncol(Qbb),nlevels(gr)) * 1:nlevels(gr))[-(nlevels(gr)+1)] +1 
          spaltEND <- c(rep(ncol(Qbb),nlevels(gr)) * 1:nlevels(gr))
          #aDD
          zeilBEG <- c(0,rep(nrow(Qbb),nlevels(gr)) * 1:nlevels(gr))[-(nlevels(gr)+1)] +1
          zeilEND <- c(rep(nrow(Qbb),nlevels(gr)) * 1:nlevels(gr))
          
          
          prmA <- c("alpha","beta","zeta","lam")
          
          
          # gehe die zeilen durch
#           galaZ=(1:nlevels(gr))[1]
#           begz=zeilBEG[1]
#           endz=zeilEND[1]
          #
          zerg <- mapply(function(begz,endz,galaZ) # speichert alle zeilen
          {
#             begsp=spaltBEG[2]
#             endsp=spaltEND[2]
            # grunr=(1:nlevels(gr))[2]
            # geht die spalten durch:
            prot <- mapply(function(begsp,endsp,grunr) # speichert je eine vollständige zeile
            {
              
              Qtemp <- Q[begz:endz,begsp:endsp] #auswahl des quadranten
              
              
              for(EACH in 1:length(design))
              {
                gala  <- design[[EACH]]
                loe   <- which(gala[galaZ,] != grunr)
                prm   <- prmA[EACH]
                
                if(length(loe) == 0)
                {
                  #Qtemp
                  next
                  
                } else {
                  
                  iaus <- paste("^.+(",paste("I",loe,collapse="|",sep=""),")",prm,sep="")
                  wo <- grep(iaus,rownames(Qtemp),value=F,perl=TRUE)
                  Qtemp[wo,] <- 0
                  
                }
              }
              Qtemp
            },begsp=spaltBEG, endsp=spaltEND, grunr=1:nlevels(gr),SIMPLIFY=FALSE)
            
            
            
          },begz=zeilBEG, endz=zeilEND,galaZ=1:nlevels(gr),SIMPLIFY=FALSE)
          
          
          zwQ <- lapply(zerg,function(XX)do.call(cbind,XX))
          Q  <- do.call(rbind,zwQ)
          onlyZ <- which(apply(Q,2,function(x)all(x == 0)))
          Q <- Q[,-onlyZ]
          
          rownames(Q) <- rwn
          colnames(Q) <- paste("eta",1:ncol(Q),sep="")
          
          #test <- Q # NEU!!!

        
        
          } else {stop("Check your input for argument: 'design' ")}
      
    }
    

    
  }
  
  return(Q)

}
