### nonparameteric test
## input is reshape object

# head(daf[[1]])
# head(my_resh$recm[[1]])


### focal group

DDF <- function(reshOBJ, wm="focal")
{
  
#   d01_r <- reshOBJ$recm[[1]]
#   daf_r <- reshOBJ$d[[1]]
#   
#   d01_f <- reshOBJ$recm[[2]]
#   daf_f <- reshOBJ$d[[2]]

# ------------------------------------------------------------------


  
faclev <- lapply(1:ncol(reshOBJ$d[[1]]), function(zi){
                    sort(unique(c(unique(reshOBJ$d[[1]][,zi]),unique(reshOBJ$d[[2]][,zi]))))
                                                  })

Ps <- mapply(function(daf,d01)
  {# groups
  corcat  <- grep("cor",colnames(d01))
  rands   <- rowSums(d01[,corcat])
  rands_f <- factor(rands,levels=0:length(corcat))
  
  groupP <- mapply(function(dd,fl)
      {#items
      mycol <- factor(daf[,dd],levels=fl)
      
      prp <- tapply(mycol,rands_f,function(ee)
              {
              prop.table(table(ee))
              })
      
      # delete NULL entries
      null_table <- table(mycol)
      null_table[] <- 0
      
      prp2 <- lapply(prp,function(inn)
              {
              if(is.null(inn)) null_table
                  else inn
              })
            },dd=1:ncol(daf),fl=faclev,SIMPLIFY=FALSE)
        
  }, daf=reshOBJ$d, d01 = reshOBJ$recm,SIMPLIFY=FALSE)
  

## und die jetzt irgendwie voneinander geschickt abziehen - ist das wunderschön oder was?


## da muss man geschickt aus der funktion vorher rausfischen




wm_all <- mapply(function(daf,d01)
{# groups
  corcat  <- grep("cor",colnames(d01))
  rands   <- rowSums(d01[,corcat])
  rands_f <- factor(rands,levels=0:length(corcat))
  
  groupP <- mapply(function(dd,fl)
  {#items
    mycol <- factor(daf[,dd],levels=fl)
    
    prp <- tapply(mycol,rands_f,function(ee)
    {
      sum(table(ee))
    })
      
    prp2 <- lapply(prp,function(inn)
    {
      if(is.na(inn)) 0
      else inn
    })
  },dd=1:ncol(daf),fl=faclev,SIMPLIFY=FALSE)
  
}, daf=reshOBJ$d, d01 = reshOBJ$recm,SIMPLIFY=FALSE)




if(wm == "focal")
{
wmc <- wm_all[[2]]  
  
}



items <- ncol(reshOBJ$d[[1]])

stdpdif <- lapply(1:items,function(alli)
    {# geht die items durch
    
      totsum <- sum(unlist(wmc[[alli]]))
    
    isws <- mapply(function(a1,a2,w)
              { # geht die randummen durch
            
              w*(a1 - a2)/totsum
              
              },a1=Ps[[1]][[alli]],a2=Ps[[2]][[alli]],w=wmc[[alli]], SIMPLIFY=TRUE)
    
    rowSums(isws)
    
    })


stdpdif
}








