reshMG <-
function(da,items=0,groups=NA,correct,design="nodif",echo=TRUE,TYPE="NRM")
{
  # d = data frame (including factors of reponses --> responses as integer beginning with 1 ending with m. contains the responses of the examinees in one of the m possible categories of item i
  #items = numerical vector which refers to the columns of the items
  # groups = numerical vector (of length = 1) which refers to the grouping variable
  # correct = vector of length(correct) == ncol(d), which contains the number of categories for each item which are counted as true
  # design = this argument refers to the function which builds the 'group design matrix'. 
  # there are 4 valid values for design:
  #1) nodif
  #2) dif1 - dif3
  
  #---------------------------
  # find items and groups ----
  #---------------------------
  
  if(is.na(groups) & any(items == 0))
  {
    d  <- da
    gr <- factor(rep("A",nrow(da)))
    items <- 1:ncol(da)
  } else if(is.na(groups)){
    d  <- da[,items]
    gr <- factor(rep("A",nrow(da)))
  } else {
    d  <- da[,items]
    gr <- da[,groups] 
  }
  
  # doing controls - its necessary - trust me
  ctrlERG <- dctrl(d,correct)


  if(ctrlERG$probl)
  {
    print(ctrlERG$problemlist)
    stop("we've got a problem")
  }
  
  # create dummy coded data
  recm <- lapply(levels(gr), function(LE)
            {
             
            dummycod <- mapply(function(well,item)
                    {
                      COL2 <- factor(ifelse(d[gr == LE,item] == well,"cor",paste("dt_",d[gr == LE,item],sep="")))
                      caa  <- model.matrix(~ -1 + COL2, data = model.frame(~ -1 + COL2, na.action=na.pass))
                      colnames(caa) <- gsub("COL2",paste("I",item,sep=""),colnames(caa))
            
                      return(caa)
                    },well=correct,item = items,SIMPLIFY = FALSE)  
            return(dummycod)
            })
  
  
  d1uc <- lapply(levels(gr), function(LE)
            {
              dummycod <- mapply(function(well,item)
              {
                COL2 <- factor(ifelse(d[gr == LE,item] == well,"cor",paste("dt_",d[gr == LE,item],sep="")))
                caa  <- model.matrix(~ -1 + COL2, data = model.frame(~ -1 + COL2, na.action=na.pass))
                colnames(caa) <- gsub("COL2",paste("I",item,sep=""),colnames(caa))
                cAA  <- ifelse(is.na(caa),0,caa)
                
                return(cAA)
              },well=correct,item = items,SIMPLIFY = FALSE)  
              return(dummycod)
            })
  

  aDD <- mapply(function(well,item)
          {
            COL2    <- factor(ifelse(d[,item] == well,"cor",paste("dt_",d[,item],sep="")))
            tabcat  <- table(COL2)
            categ   <- levels(COL2)
            anz_cat <- length(categ)
            addit   <- list(tabcat=tabcat,categ=categ,anz_cat=anz_cat)
            #
            return(addit)
          },well=correct,item = items,SIMPLIFY = FALSE)


  
  #recm <- lapply(dummycod,function(x)x[[1]])
  #aDD  <- lapply(dummycod,function(x)x[[2]])
  #d1uc <- lapply(recm,function(x) lapply(x, function(xx1) ifelse(is.na(xx1),0,xx1)))
  
  cat("data = reshaped\n")
  
  # descriptives
  coluN <- max(sapply(d,function(x)length(table(x))))
  relev <- data.frame(lapply(d,function(x) factor(x,levels=1:coluN,labels=1:coluN)))
  
  
  
  if(is.na(groups))
  {
    absF1 <- sapply(relev,function(TA)table(TA,useNA="always"))
    rownames(absF1)[1:coluN] <- paste("category",1:coluN,sep="")  
  } else {
    absF1 <- lapply(relev,function(TA)
    {
      tt1 <- table(TA,gr,useNA="always")
      rownames(tt1)[1:coluN] <- paste("category",1:coluN,sep="")
      tt1
    })
  }
 
  # controlling the input design !!
  if(is.list(design))
      {
      design <- ctrl_design(design=design,aDD=aDD,gr=gr,TYPE=TYPE)  
      }
      

  # create Q matrix

      gdema <- grDM(aDD,gr,design,TYPE=TYPE)
 
  
  if(echo){print(absF1)}
  
  cat("group information added \n")
  reshret <- list(recm=recm,aDD=aDD,d=d,gr=gr,Qmat=gdema,d1uc=d1uc,design=design)
  
  # classify
  if(TYPE=="NRM")
      {
        class(reshret) <- "reshNRM" 
      } else if(TYPE=="NLM") {
                              class(reshret) <- "reshNLM"  
                             }
      
  return(reshret)
}
