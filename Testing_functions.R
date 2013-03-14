# now lets try different scenarios
####----------------------------------------------------------------


### EXECUTE THE FUNCTIONS OF NLM
# --------------------------------------------------

# Make a Parlist for the case of 4 Items

Item1 <- c(1,0,c(-0.5,0.3,0.2),c(-0.5,-0.3,0.8))
#Item1 <- c(1,-2,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
#Item1 <- c(1,-2,c(-0.5,-0.2,0.4,0.3),c(-1,0.3,0.1,0.6))
names(Item1) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))
#names(Item1) <- c("a","b",paste("zeta1",1:4,sep=""),paste("lamb",1:4,sep=""))

Item2 <- c(1,-0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item2) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item3 <- c(1,0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item3) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item4 <- c(1,0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item4) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

ParList <- list(Item1=Item1,Item2=Item2,Item3=Item3,Item4=Item4)


set.seed(11)
pp <- rnorm(1000)
#pp <- rnorm(10000)
#set.seed(12)
sim.nlm.1 <- NLM.sim(ParList,pp)

head(sim.nlm.1)


sim.nlm.1[c(2,15,500:509,799:803),1] <- NA
sim.nlm.1[c(12,19,20,900:940),4] <- NA


# reshape
dre <- resh(sim.nlm.1,c(1,1,1,1))
d1 <- dre[[1]]



# starting values
startV <- startV_nlm(dre)

ulstv <- startV$ulstv
SKEL <- startV$stwm1

# NODES
library(statmod)
quads <- gauss.quad.prob(15,"normal",mu=0)


#!
erg_nlm1 <- NestLM(dre)

erg_nlm1$parest

#erg_estep <- Enlm(Ulstv,d1,quads,SKEL)
#names(erg_estep)



###########################################
## missing values #########################
###########################################


#  ++++ EM Algorithm
# ======================================
#d1    <- dre[[1]]
OLD   <- 0
ZAEHL <- 1
PARS  <- ulstv
system.time(
  repeat
  {
    cat(ZAEHL,"\r") 
    
    #E ****
    erg_estep <- Enlm(PARS,d1,quads,SKEL)
    
    ##M ****
    oerg <- optim(par=PARS,fn=ZFnlm,gr=de1nlm,erg_estep=erg_estep,quads=quads,
                  SKEL=SKEL,control=list(fnscale=-1,maxit=70),method="BFGS",hessian=F)
    
    #oerg <- optim(par=PARS,fn=ZFnlm,erg_estep=erg_estep,quads=quads,
    #              SKEL=SKEL,control=list(fnscale=-1,maxit=50),method="BFGS",hessian=TRUE)
    
    if(abs(OLD - abs(oerg$value)) <= 0.0001 | ZAEHL > 1000){break}
    
    OLD <- abs(oerg$value)
    PARS <- oerg$par
    ZAEHL <- ZAEHL + 1 
  })
# es wird berechnet --> die ergebnisse unterscheiden sich (logischerweise) leicht von den anderen.
# habe am 08.02.2012 nochmal den e-step einzeln ausgeführt um zu überprüfen inwiefern mit NAs korrekt ungegangen wird, und es ist alles plausibel und dürfte passen!







set.seed(11)
pp <- rnorm(1000)
#pp <- rnorm(10000)
set.seed(12)
sim.nlm.1 <- NLM.sim(ParList,pp)

head(sim.nlm.1)


#sim.nlm.1[c(2,15,500:509,799:803),1] <- NA
#sim.nlm.1[c(12,19,20,900:940),4] <- NA


# reshape
dre <- resh(sim.nlm.1,c(1,1,1,1))


#NestLM <- function(reshOBJ,METH="BFGS",exac=0.001,maxitEM=500,nQUAD=15,...)

#!
erg_nlm1 <- NestLM(dre)

erg_nlm1$parest

plot(density(erg_nlm1$eap_est,bw=1))

dctrl <- data.frame(sim.nlm.1,eap=erg_nlm1$eap_est)


dctrl[apply(sim.nlm.1,1,function(x) sum(x==1) == 4),]
dctrl[apply(sim.nlm.1,1,function(x) sum(x==2) == 4),]

dctrl[which.max(erg_nlm1$eap_est),]
dctrl[which.min(erg_nlm1$eap_est),]



dctrl2 <- data.frame(sim.nlm.1,eap=thsc)

dctrl2[apply(sim.nlm.1,1,function(x) sum(x==2) == 4),]


dctrl2[which.min(dctrl2$eap),]
dctrl2[which.max(dctrl2$eap),]


plot(erg_nlm1$eap_est,ppo$par,ylim=c(-10,10))



summary(erg_nlm1$eap_est)

nlmOBJ <- erg_nlm1

nlmplot1 <- plot(erg_nlm1,1:2)
nlmplot1






########################################
### PP experiment
########################################




#THET <- rep(0.1,1000)

# MLE person parameter estimation
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PPnlm <- function(THET,Ulstv,d1,SKEL)
{
  relstv <- relist(Ulstv,SKEL)
  
  
  nrme1 <- mapply(function(pitem,daten) # rattert pro item durch
  {
    # 2PL Part
    # ----------------------
    Km      <- matrix(c(rep(1,length(THET)),THET),ncol=2)
    
    abpar <- c(-pitem[2],-pitem[1]) # hat sich gegenüber der vorgängerversion geändert - ist jetzt so wie bei yss beschrieben -- muss wohl so sein, um die selbe polung zu haben wie NRM
    solit <- 1/ (1 + exp(Km %*% abpar))
    
    dosolit <- 1-solit
    tplpart <- cbind(solit,dosolit)
    
    # NRM Part
    # ----------------------
    pitemNRM <- pitem[-(1:2)]
    
    ZQstern <- coP_nrm(pitemNRM,Km)
    
    ZQstern_nlm <- ZQstern * as.vector(dosolit)
    
    # now the data
    # ---------------------
    
    #corrcat <- daten[,1]
    #distr   <- t(daten[,-1])
    ZQstern_all <- cbind(solit,ZQstern_nlm)
    
    ZQdstern <- rowSums(ZQstern_all * daten)
    # missing values remain as NA in the Data-Structure
    
    sum(log(ZQdstern))
  },pitem=relstv,daten=d1,SIMPLIFY = TRUE)
  
  
  
  sum(nrme1)
}

THET <- rep(0.1,1000)
ppo <- optim(par=THET,fn=PPnlm,Ulstv=erg_nlm1$optim_obj$par,
             d1=d1,SKEL=SKEL,control=list(fnscale=-1,maxit=100),method="BFGS",hessian=T)


ppo$par
ppdf <- data.frame(sim.nlm.1,ppmle=ppo$par)

ppdf[which.min(ppo$par),]



vgl <- data.frame(dctrl,ppmle=ppo$par,zaehler=test)

vgl[which.min(ppo$par),]






##################################################################################################
############################################################################################
##### recreating the functions to multigroup and Q Matrix things
############################################################################################
##################################################################################################




Item1 <- c(1,-2,c(-0.5,0.3,0.2),c(-0.5,-0.3,0.8))
names(Item1) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item2 <- c(1,-1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item2) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item3 <- c(1,0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item3) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item4 <- c(1,1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item4) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item5 <- c(1,2,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item5) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

ParList <- list(Item1=Item1,Item2=Item2,Item3=Item3,Item4=Item4,Item5=Item5)


set.seed(555555)
pp <- rnorm(1000)
#pp <- rnorm(10000)
#set.seed(12)
sim.nlm.1 <- NLM.sim(ParList,pp)

head(sim.nlm.1)

# reshape
#dre <- resh(sim.nlm.1,c(1,1,1,1,1))

# now use the new reshMG function which is identically adopted from the NRM folder and copied to the shared functions folder. so lets whether this gives useful results for further processing.

dre_new <- reshMG(sim.nlm.1,items=1:5,groups=NA,correct=c(1,1,1,1,1), TY="NLM")

# first, we have to checkout the control function, and check wheater this controls are compatible to the NLM model

#correct <- c(1,1,1,1)
#d <- sim.nlm.1

###### change the Qmatrix function to fit in the logic of NLM


# 1: ändern der zeilennamen in Q
# 2: ausprobieren ob auch dasselbe rauskommt - d.h. umbauen der logLIK etc. ... für die neue Q matrix
# 3: wenn da alles rennt dann commit und weitere umabauten bei Qmatrix und design etc. 


reshOBJ <- dre_new
# astart=1
# bstart=0.5
# gamstart=0.1
# xistart=-0.1

# check out how the nrm starting values look like

startOBJ <- startV_nlmMG(reshOBJ)

#### jump into Estep

Ulstv <- startOBJ$ulstv

quads <- quadIT()
# ws end XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxx





esteperg <- Enlm(Ulstv=Ulstv,reshOBJ=reshOBJ,startOBJ=startOBJ,quads=quads)
# Result looks like:
## $group1
### $riq_querG
### $riqv_querG
# ETC.
erg_estep <- esteperg


zferg <- ZFnlm(Ulstv=Ulstv,erg_estep=erg_estep,startOBJ=startOBJ,reshOBJ=reshOBJ,quads=quads)
# scheint zu passen

# ersten Ableitungen
#de1nlm <- function(Ulstv,erg_estep,startOBJ,reshOBJ,quads)
ersteABL_newGfunc <- de1nlm(Ulstv=Ulstv,erg_estep=erg_estep,startOBJ=startOBJ,reshOBJ=reshOBJ,quads=quads)
# it seems great

mueNLM(Ulstv=Ulstv,startOBJ=startOBJ,reshOBJ=reshOBJ,quads=quads,sigmaest=FALSE)

### change the starting values to reasonable ones:

test <- startV_nlmMG(reshOBJ)


## also problem => die ersten ableitungen

### einmal die alten ersten ableitungen umbauen mit der neuen Q matrix um zu probieren ob das zu denselben gschichten führt!



###>>>>> jetzt mal EM-alg durchlaufen lassen (in mainfunc) und schauen was rauskommt, und mal heuristisch überprüfen:
#############################################################################
test <- as.vector(reshOBJ$Qmat %*% ESTlist[[1]])

test <- as.vector(reshOBJ$Qmat %*% oerg$par)

test[2:4] - mean(test[2:4])
test[6:8] - mean(test[6:8])

test[10:12] - mean(test[10:12])
test[14:16] - mean(test[14:16])

### interessant. es scheint noch ein fehler in den ersten Ableitungen zu sein.


#   Item1 <- c(1,-2,c(-0.5,0.3,0.2),c(-0.5,-0.3,0.8))
#   names(Item1) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))
#   
#   Item2 <- c(1,-1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
#   names(Item2) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))
#   
#   Item3 <- c(1,0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
#   names(Item3) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))
#   
#   Item4 <- c(1,1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
#   names(Item4) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))
#   
#   Item5 <- c(1,2,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
#   names(Item5) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))
#   


# rewriting Q matrix for NLM
######################################### ######################################

# dre_new from line 287

aDD <- dre_new$aDD
gr <- dre_new$gr
TYPE <- "NLM"


### but now we need at least 2 groups to build the new Q function


perp1 <- rnorm(1000,0,1)
perp2 <- rnorm(1000,1,1)

simdat1 <- NLM.sim(ParList,perp1)
simdat2 <- NLM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:1000,simdat1)
simdat2 <- data.frame(ID=1001:2000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=T)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=1000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)


items=1:5
correct=rep(1,5)
groups=6
design = "nodif"
echo=TRUE
da <- DAT1

#reshmulti1 <- reshMG(DAT1,,)

# design for NLM
# now we create a list which contains 4 matrices to define the dif design for the NLM model
# 1 --> beta pars
# 2 --> alpha pars
# 3 --> zetapars
# 4 --> lambdapars

design1 <- designTemp(2,5,"NLM")
design <- design1




nlmnodif <- grDM(aDD=aDD,gr=gr,design="nodif",TYPE="NLM")
write.csv(nlmnodif,file="/home/manuel/Dokumente/Manuel/packages/mcIRT/NLM/nlm_Funktionen_Tests/Qfunctest/Qnlm.csv")


design[[1]][2,3:5] <- 2

nlmdif1 <- grDM(aDD=aDD,gr=gr,design=design,TYPE="NLM")
write.csv(nlmdif1,file="/home/manuel/Dokumente/Manuel/packages/mcIRT/NLM/nlm_Funktionen_Tests/Qfunctest/Qnlmdif1.csv")



design[[1]][2,3:4] <- 2
design[[2]][2,c(1,3,5)] <- 2
design[[3]][2,c(2:5)] <- 2

nlmdif2 <- grDM(aDD=aDD,gr=gr,design=design,TYPE="NLM")
write.csv(nlmdif2,file="/home/manuel/Dokumente/Manuel/packages/mcIRT/NLM/nlm_Funktionen_Tests/Qfunctest/Qnlmdif2.csv")

### control thing 4 design


design[[1]][2,3] <- 3
design <- design[-3]


ctrl_design(design=design,aDD=aDD,gr=gr,TYPE=TYPE)



################################### now testing the multigroup estimation:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



perp1 <- rnorm(1000,0,1)
perp2 <- rnorm(1000,1,1)

simdat1 <- NLM.sim(ParList,perp1)
simdat2 <- NLM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:1000,simdat1)
simdat2 <- data.frame(ID=1001:2000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=T)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=1000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)

reshOBJ <- reshMG(DAT1,items=1:5,groups=6,correct=rep(1,5),design="nodif",echo=TRUE,TYPE="NLM")


############## testing different parametrizations
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+


Item1 <- c(1,-2,c(-0.5,0.3,0.2),c(-0.5,-0.3,0.8))
names(Item1) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item2 <- c(1,-1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item2) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item3 <- c(1,0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item3) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item4 <- c(1,1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item4) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item5 <- c(1,2,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item5) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

ParList <- list(Item1=Item1,Item2=Item2,Item3=Item3,Item4=Item4,Item5=Item5)



set.seed(741258)
perp1 <- rnorm(1000,0,1)

simdat1 <- NLM.sim(ParList,perp1)

reshOBJ <- reshMG(simdat1,items=1:5,groups=NA,correct=rep(1,5),design="nodif",echo=TRUE,TYPE="NLM")

test <- reshOBJ$Qmat

test[grep("(alpha|beta)",rownames(test)),] <- test[grep("(alpha|beta)",rownames(test)),]*(-1)


reshOBJ$Qmat <- test

erg.allMIN1 <- nlm(reshOBJ=reshOBJ)


etaerg <- as.vector(reshOBJ$Qmat %*%  erg.allMIN1$etapar)


etaerg[grep("(beta)",rownames(test))] - mean(etaerg[grep("(beta)",rownames(test))])
# betas sind verdreht
etaerg[grep("(alpha)",rownames(test))]
# alphas sind korrekt

etaerg[2:4] - mean(etaerg[2:4])
# zetas sind korrekt

etaerg[6:8] - mean(etaerg[6:8])
# lambdas sind verdreht



beta1 <- 0.2
alpha1 <- 1.2
theta1 <- -0.7


1/(1+exp(beta1 + alpha1*theta1))

1/(1+exp(-beta1 - alpha1*theta1))

exp(beta1 + alpha1*theta1)/(1+exp(beta1 + alpha1*theta1))
exp(-beta1 - alpha1*theta1)/(1+exp(-beta1 - alpha1*theta1))


# now changing the starting values to zero for each parameter
#-------------------------------------------------------------------------

erg.allMIN2 <- nlm(reshOBJ=reshOBJ,etastart=0)

etaerg <- as.vector(reshOBJ$Qmat %*%  erg.allMIN2$etapar)


etaerg[grep("(beta)",rownames(test))] - mean(etaerg[grep("(beta)",rownames(test))])
# betas sind verdreht
etaerg[grep("(alpha)",rownames(test))]
# alphas sind verdreht

etaerg[2:4] - mean(etaerg[2:4])
# zetas sind korrekt
etaerg[6:8] - mean(etaerg[6:8])
# lambdas sind korrekt


# introducing the new nlm_twoplpart inkl zero starting values and Q-matrix is -1 everywhere

erg.allMIN3 <- nlm(reshOBJ=reshOBJ,etastart=0)


etaerg <- as.vector(reshOBJ$Qmat %*%  erg.allMIN3$etapar)


etaerg[grep("(beta)",rownames(test))] - mean(etaerg[grep("(beta)",rownames(test))])
# betas passen
etaerg[grep("(alpha)",rownames(test))]
# alphas sind verdreht

etaerg[2:4] - mean(etaerg[2:4])
# zetas sind korrekt
etaerg[6:8] - mean(etaerg[6:8])
# lambdas sind verdreht


############### Qmatrix auf 1 überall - rest bleibt gleich

test <- reshOBJ$Qmat

#test[grep("(alpha|beta)",rownames(test)),] <- test[grep("(alpha|beta)",rownames(test)),]*(-1)


reshOBJ$Qmat <- test * (-1)

erg.allMIN4 <- nlm(reshOBJ=reshOBJ)


etaerg <- as.vector(reshOBJ$Qmat %*%  erg.allMIN4$etapar)


etaerg[grep("(beta)",rownames(test))] - mean(etaerg[grep("(beta)",rownames(test))])
# betas passen
etaerg[grep("(alpha)",rownames(test))]
# alphas passen

etaerg[2:4] - mean(etaerg[2:4])
# zetas passen
etaerg[6:8] - mean(etaerg[6:8])
# lambdas passen


############################################ start = 0

erg.allMIN5 <- nlm(reshOBJ=reshOBJ,etastart=0)

etaerg <- as.vector(reshOBJ$Qmat %*%  erg.allMIN5$etapar)


etaerg[grep("(beta)",rownames(test))] - mean(etaerg[grep("(beta)",rownames(test))])
# betas passen
etaerg[grep("(alpha)",rownames(test))]
# alphas verdreht

etaerg[2:4] - mean(etaerg[2:4])
# zetas passen
etaerg[6:8] - mean(etaerg[6:8])
# lambdas verdreht



############################################ start = 0.1

erg.allMIN6 <- nlm(reshOBJ=reshOBJ,etastart=0.1)

etaerg <- as.vector(reshOBJ$Qmat %*%  erg.allMIN6$etapar)


etaerg[grep("(beta)",rownames(test))] - mean(etaerg[grep("(beta)",rownames(test))])
# betas passen
etaerg[grep("(alpha)",rownames(test))]
# alphas passen

etaerg[2:4] - mean(etaerg[2:4])
# zetas passen
etaerg[6:8] - mean(etaerg[6:8])
# lambdas verdreht

  erg.allMIN6$Parcent$nrmpar

##################################################
# rewritten Q matrix

reshOBJ2 <- reshMG(simdat1,items=1:5,groups=NA,correct=rep(1,5),design="nodif",echo=TRUE,TYPE="NLM")

reshOBJ2$Qmat

system.time(erg.allMIN7 <- nlm(reshOBJ=reshOBJ2))
erg.allMIN7$Parcent

################################################################
####################################################################################



direc1 = "/media/82190a37-7663-4e25-9ced-e9e55a6a00e6/homelink/Manuel/packages/mcIRT/NLM"
direc2 = "/media/82190a37-7663-4e25-9ced-e9e55a6a00e6/homelink/Manuel/packages/mcIRT/shared_functions"
direc <- list(direc1,direc2)

loadALL(direc)



Item1 <- c(1,-2,c(-0.5,0.3,0.2),c(-0.5,-0.3,0.8))
names(Item1) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item2 <- c(1,-1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item2) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item3 <- c(1,0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item3) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item4 <- c(1,1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item4) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item5 <- c(1,2,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item5) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

ParList <- list(Item1=Item1,Item2=Item2,Item3=Item3,Item4=Item4,Item5=Item5)


names(ParList) <- paste("item",1:5,sep="")
# Add names to each list element, because this names are added to the columns in the simulated data.frame.



perp1 <- rnorm(1000,0,1)
perp2 <- rnorm(1000,1,1)

simdat1 <- NLM.sim(ParList,perp1)
simdat2 <- NLM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:1000,simdat1)
simdat2 <- data.frame(ID=1001:2000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=T)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=1000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)

des1 <- designTemp(ngru=2,nit=5,TYPE="NLM")

reshdat <- reshMG(DAT1,items=1:5,groups=6,correct=rep(1,5),TYPE="NLM",echo=FALSE,design="nodif")
ergnlm1 <- nlm(reshdat)

ergnlm1$Parcent

########################################################
########################################################################

direc1 = "/home/manuel/Dokumente/Manuel/packages/mcIRT/NLM"
direc2 = "/home/manuel/Dokumente/Manuel/packages/mcIRT/shared_functions"
direc <- list(direc1,direc2)
loadALL(direc=direc)


RESnlm <- ergnlm1
# names(RESnlm)[8] <- "ZLpar"
# names(RESnlm$erg_distr) <- c("mean_est","sig_est")
# 
# reshOBJ <- reshdat
# startOBJ <- startV_nlmMG(reshOBJ=reshOBJ,etastart="aut")
# quads <- quadIT(ngr=2)
# ergmue <- mueNLM(Ulstv=startOBJ$ulstv,reshOBJ,startOBJ,quads,T,T)
# ergmue$errmat



summary(RESnlm)
plot(RESnlm)
print(RESnlm)

des1 <- designTemp(2,5,TYPE="NLM")

des1[[2]][2,1:3] <- 2 


reshdat <- reshMG(DAT1,items=1:5,groups=6,correct=rep(1,5),TYPE="NLM",echo=FALSE,design=des1)
ergnlm1 <- nlm(reshdat)

##### build package:

direc1 = "/home/manuel/Dokumente/Manuel/packages/mcIRT_preRelease/NLM"
direc2 = "/home/manuel/Dokumente/Manuel/packages/mcIRT_preRelease/shared_functions"
direc3 = "/home/manuel/Dokumente/Manuel/packages/mcIRT_preRelease/NRM"
direc <- list(direc1,direc2,direc3)
loadALL(direc=direc)

rm(list=c("direc1","direc2","direc3","direc"))

package.skeleton(name="mcIRT",path="/mnt/datasdisk/homelink/Manuel/packages")




#################################################################




Item1 <- c(1,-2,c(-0.5,0.3,0.2),c(-0.5,-0.3,0.8))
names(Item1) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item2 <- c(1,-1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item2) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item3 <- c(1,0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item3) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item4 <- c(1,1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item4) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item5 <- c(1,2,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item5) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

ParList <- list(Item1=Item1,Item2=Item2,Item3=Item3,Item4=Item4,Item5=Item5)


names(ParList) <- paste("item",1:5,sep="")
# Add names to each list element, because this names are added to the columns in the simulated data.frame.



perp1 <- rnorm(1000,0,1)
perp2 <- rnorm(1000,1,1)

simdat1 <- NLM.sim(ParList,perp1)
simdat2 <- NLM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:1000,simdat1)
simdat2 <- data.frame(ID=1001:2000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=T)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=1000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)

des1 <- designTemp(ngru=2,nit=5,TYPE="NLM")

reshdat <- reshMG(DAT1,items=1:5,groups=6,correct=rep(1,5),TYPE="NLM",echo=FALSE,design="nodif")
ergnlm1 <- nelm(reshdat)

plot(ergnlm1,fromto=c(-6,6))



ergnlm1$last_mstep$hessian


TIMEme <- nelm(reshdat,EMmax=8)

######################################## NRM

library(mcIRT)


NUMBI <- 5

ParList <- lapply(1:NUMBI,function(x)
{
  Item1 <- c(c(-2,-1,1,2),c(-1.2,0.3,0.2,0.7))
  names(Item1) <- c(paste("zeta",1:4,sep=""),paste("lamb",1:4,sep=""))
  Item1
})


names(ParList) <- paste("item",1:NUMBI,sep="")
# Add names to each list element, because this names are added to the columns in the simulated data.frame.



perp1 <- rnorm(1000,0,1)
perp2 <- rnorm(1000,1,1)

simdat1 <- NRM.sim(ParList,perp1)
simdat2 <- NRM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:1000,simdat1)
simdat2 <- data.frame(ID=1001:2000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=T)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=1000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)

# hier wird alles mit den 2 gruppen als "nodif" geschätzt
reshdat <- reshMG(DAT1,items=1:5,groups=6,correct=rep(1,5))

system.time(ergnrmN <- nrm(reshdat,EMmax=200,sigmaest=TRUE))

plot(ergnrmN,1000,c(-8,8))



#########################



erg.allMIN1 <- nlm(reshOBJ=reshOBJ)


###########################################################################################
############################################## creating the second derivates for nrm
###################################################################################


NUMBI <- 5

ParList <- lapply(1:NUMBI,function(x)
{
  Item1 <- c(c(-2,-1,1,2),c(-1.2,0.3,0.2,0.7))
  names(Item1) <- c(paste("zeta",1:4,sep=""),paste("lamb",1:4,sep=""))
  Item1
})
names(ParList) <- paste("item",1:NUMBI,sep="")
perp1 <- rnorm(1000,0,1)
simdat1 <- NRM.sim(ParList,perp1)


reshdat <- reshMG(simdat1,items=1:5,correct=rep(1,5))




# starting object:
erg_start <- mcIRT:::startV_nrmMG(reshdat)
QUADS <- mcIRT:::quadIT()
PREVinp <- NA
ergEstep <- mcIRT:::Enrm(erg_start$ulstv, reshdat, erg_start, QUADS, PREVinp)


#NRM_2deriv <- function(Ulstv,riqv_quer,quads,SKEL)

Ulstv <- erg_start$ulstv
riqv_quer <- ergEstep
reshOBJ  <- reshdat
quads <- QUADS  
startOBJ <- erg_start   
  
################################

testerg <- nrm(reshdat)
Ulstv <- testerg$last_mstep$par
riqv_quer <- testerg$last_estep
  
testerg$last_mstep$hessian  
  
  #####################



perp1 <- rnorm(10000,0,1)
simdat1 <- NRM.sim(ParList,perp1)


reshdat <- reshMG(simdat1,items=1:5,correct=rep(1,5))

ti1 <- system.time(testerg1 <- nrm(reshdat, dooptim=TRUE, EMmax=1000))
ti2 <- system.time(testerg2 <- nrm(reshdat, dooptim=FALSE, EMmax=1000))

testerg1$last_mstep$value
testerg2$last_mstep$value


############### UMBAU VON nrm #########
#######################################

ctrl <- list()
ctrl$nodes <- 22
ctrl$verbose <- FALSE

ctrl <- list()
ctrl$verbose <- FALSE
ctrl$nodes <- 22
ctrl$wurstbrot <- "wurst"


########################################################




NUMBI <- 5

ParList <- lapply(1:NUMBI,function(x)
{
  Item1 <- c(c(-2,-1,1,2),c(-1.2,0.3,0.2,0.7))
  names(Item1) <- c(paste("zeta",1:4,sep=""),paste("lamb",1:4,sep=""))
  Item1
})
names(ParList) <- paste("item",1:NUMBI,sep="")
perp1 <- rnorm(1000,0,1)
simdat1 <- NRM.sim(ParList,perp1)


reshdat <- reshMG(simdat1,items=1:5,correct=rep(1,5))


changedNRM1a <- nrm(reshdat)
changedNRM1 <- nrm(reshdat, ctrl=list(nodes=8, EMmax=400, eu1=22, xml=55))

# check out the information functions
ESTlist <- changedNRM1

changedNRM1b <- nrm(reshdat, ctrl=list(nodes=8, EMmax=400, NRexac=0.0001))


mcIRT:::plotInrm(changedNRM1)


###################################### NLM second deriv und laufzeitverhalten




Item1 <- c(1,-2,c(-0.5,0.3,0.2),c(-0.5,-0.3,0.8))
names(Item1) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item2 <- c(1,-1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item2) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item3 <- c(1,0,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item3) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item4 <- c(1,1,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item4) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

Item5 <- c(1,2,c(-0.5,-0.3,0.8),c(-0.5,0.3,0.2))
names(Item5) <- c("a","b",paste("zeta1",1:3,sep=""),paste("lamb",1:3,sep=""))

ParList <- list(Item1=Item1,Item2=Item2,Item3=Item3,Item4=Item4,Item5=Item5)


names(ParList) <- paste("item",1:5,sep="")
# Add names to each list element, because this names are added to the columns in the simulated data.frame.



perp1 <- rnorm(5000,0,1)
perp2 <- rnorm(5000,1,1)

simdat1 <- NLM.sim(ParList,perp1)
simdat2 <- NLM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:5000,simdat1)
simdat2 <- data.frame(ID=5001:10000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=T)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=5000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)

des1 <- designTemp(ngru=2,nit=5,TYPE="NLM")

reshdat <- reshMG(DAT1,items=1:5,groups=6,correct=rep(1,5),TYPE="NLM",echo=FALSE,design="nodif")
ergnlm1 <- nelm(reshdat)

plot(ergnlm1,fromto=c(-6,6))



ergnlm1$last_mstep$hessian


TIMEme <- nelm(reshdat,EMmax=8, sigmaest=TRUE)

# Wenn die zweiten Ableitungen drinnenstehen, dann ist
# der ESTEP das langsamste (bzw. eben die berechnung von MUE) !!
# --> das w?re es wert in Rcpp auszulagern

system.time(TIMEme2 <- nelm(reshdat,ctrl=list(EMmax=500)))
#system.time(TIMEme2a <- nelm(reshdat,ctrl=list(EMmax=500,dooptim=TRUE)))



#####################################################################################

plotINF(TIMEme2)



NUMBI <- 5

ParList <- lapply(1:NUMBI,function(x)
{
  Item1 <- c(c(-2,-1,1,2),c(-1.2,0.3,0.2,0.7))
  names(Item1) <- c(paste("zeta",1:4,sep=""),paste("lamb",1:4,sep=""))
  Item1
})


names(ParList) <- paste("item",1:NUMBI,sep="")
# Add names to each list element, because this names are added to the columns in the simulated data.frame.



perp1 <- rnorm(5000,0,1)
perp2 <- rnorm(5000,1,1)

simdat1 <- NRM.sim(ParList,perp1)
simdat2 <- NRM.sim(ParList,perp2)

simdat1 <- data.frame(ID=1:5000,simdat1)
simdat2 <- data.frame(ID=5001:10000,simdat2)

simdatalla <- merge(simdat1,simdat2,all=T)
simdatall  <- simdatalla[,-1]

head(simdatall)
gruAB <- factor(rep(c("A","B"),each=5000))

DAT1 <- data.frame(simdatall,ABgroup = gruAB)

head(DAT1)

# hier wird alles mit den 2 gruppen als "nodif" geschätzt
reshOBJ <- reshMG(DAT1,items=1:NUMBI,groups=NUMBI+1,correct=rep(1,NUMBI))


########## testen im speedup branch

test <- nrm(reshOBJ)


test2 <- nrm(reshOBJ)






