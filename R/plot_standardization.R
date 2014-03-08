plot.standardization <- function(x, type="rmwsd",...)

{
  
  
  if(type == "rmwsd")
  {
    par("ask"=TRUE)

    
    for(pls in 1:length(x$varpp))
    {
    
      XYL <- max(sqrt(x$varpp[[pls]]),abs(x$stdpdif[[pls]]))

      if(XYL < 0.3)
        {
        XYL <- 0.3  
        }
      
    plot(x$stdpdif[[pls]],sqrt(x$varpp[[pls]]),xlim=c(-XYL,XYL),ylim=c(0,XYL),xaxs="i",yaxs="i",pch=4,col="red",xlab="standardized difference",ylab="SD of difference", main=names(x$stdpdif)[pls],...)
    #symbols(rep(0,6),rep(0,6),circles=c(0.05,0.075,0.1,0.15,0.2,0.3),add=TRUE, inches=FALSE,lty=2)
    symbols(rep(0,3),rep(0,3),circles=c(0.08,0.16,0.24),add=TRUE, inches=FALSE,lty=2)
    text(x$stdpdif[[pls]],sqrt(x$varpp[[pls]]),names(x$stdpdif[[pls]]),pos=4, offset=0.3)
    points(x$stdpdif[[pls]],sqrt(x$varpp[[pls]]),cex=1.4,col="red")
    abline(v=seq(-1,1,0.1),lty=2)
    
    }
  
  } else if(type == "catwise")
      {
      
     par("ask"=TRUE)
      for(it in 1:length(x$stdpdif)) # items
      {


        whereobs <- apply(x$Ps[[1]][[it]],2,function(w) !all(log(w) < -30))
        whereobs2 <- apply(x$Ps[[2]][[it]],2,function(w) !all(log(w) < -30))
        xen   <- as.numeric(as.character(unique(c(names(whereobs)[whereobs],names(whereobs2)[whereobs2]))))
        whichone <- sum(whereobs) >= sum(whereobs2)
        
        grada <- x$Ps[[1]][[it]]
        
        grada[,!whereobs] <- NA
        if(whichone)
          {
          grada <- grada[,whereobs]  
          } else 
            {
              grada <- grada[,whereobs2] 
            }
        
        # lines for group 1
        plot(xen,grada[1,],type="l",ylim=c(0,1),xlab="score",ylab="prob", main=paste("Item",it))
        text(xen,grada[1,],labels=rownames(gradaf)[1])
        for(i in 2:nrow(grada))
          {
          text(xen,grada[i,],labels=rownames(gradaf)[i])
          lines(xen,grada[i,])
          }

        
        
        grada <- x$Ps[[2]][[it]]
        
        grada[,!whereobs2] <- NA
        if(whichone)
        {
          grada <- grada[,whereobs]  
        } else 
        {
          grada <- grada[,whereobs2] 
        }
        
        # lines for group 2
        for(i in 1:nrow(grada))
        {
          text(xen,grada[i,],labels=rownames(gradaf)[i],col="red")
          lines(xen,grada[i,],col="red",lty=2)
        }
        
        
        
        
      }
    
        
      }
  
  
  
  

  
}
