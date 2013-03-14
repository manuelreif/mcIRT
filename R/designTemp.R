designTemp <-
function(ngru,nit,TYPE="NLM")
# this is a function to create a template for the "design" argument in grDM (which is used by reshMG)
  {
  if(TYPE=="NLM")
      {
        designL <- lapply(1:4,function(x)
                    {
                    matrix(1,nrow=ngru,ncol=nit)  
                    })
        
        names(designL) <- c("alpha","beta","zeta","lambda")
        
      }
  
  if(TYPE=="NRM")
      {
        designL <- lapply(1:2,function(x)
                    {
                      matrix(1,nrow=ngru,ncol=nit)  
                    })
        
        names(designL) <- c("zeta","lambda")
        
      }
  
  
  designL  
}
