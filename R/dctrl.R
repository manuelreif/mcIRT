dctrl <-
function(d,correct)
{
  problemlist <- list()  
  a <- 1  
  probl <- FALSE
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ dataframe
  if(!is.data.frame(d))
  {
    problemlist[[a]] <- "data object is not a data.frame!"
    a <- a + 1
    probl <- TRUE
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ minimum
  # minimum <- sapply(d,function(x)min(as.numeric(x),na.rm=T) == 1)
  minimum <- sapply(d,function(x)min(as.numeric(x),na.rm=T) == 0)
  if(!any(minimum))
  {
    problemlist[[a]] <- "minimum != 1 in at least one column!"
    a <- a + 1
    probl <- TRUE
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ number of cat
  
  un <- sapply(d,function(x)length(table(x))> 1)
  if(!any(un))
  {
    problemlist[[a]] <- "there are variables containing only 1s!"
    a <- a + 1
    probl <- TRUE
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ correct_1
  
  leco <- length(correct) == ncol(d)
  if(!leco)
  {
    problemlist[[a]] <- "length(correct) != ncol(d)"
    a <- a + 1
    probl <- TRUE
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ correct_2
  
  leco2 <- mapply(function(A,B) B %in% A ,A=d,B=correct)
  if(!any(leco2))
  {
    problemlist[[a]] <- "the category of correct response was not observed (in one or more variables)!"
    probl <- TRUE
  }
  
  list(problemlist=problemlist,probl=probl)
}
