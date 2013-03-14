quadIT <-
function(nodes=14, mu=0, sigma=1,absrange=5, ngr=1) 
{
  nwpgru <- mapply(function(gruu,mug,si)
  {
    
    if(length(nodes)==1)
    {
      
      if(nodes < 7) 
      {
        nodes <- 7
        cat("#nodes == too small --> set to 7\n")
        absrange <- 4
      }
      
      quadP <- seq(absrange*(-1), absrange, length.out = nodes) 
      
      quadweight  <- dnorm(quadP)
      quadweight1 <- quadweight/sum(quadweight)
      quadP_shift <- quadP*si + mug
      
    } else  {
      quadP <- nodes - mean(nodes)
      
      quadweight  <-  dnorm(quadP)
      quadweight1 <- quadweight/sum(quadweight)
      quadP_shift <- quadP*si + mug
    }
    
    list(nodes=quadP_shift,weights=quadweight1)
  },gruu=1:ngr,mug=mu,si=sigma,SIMPLIFY=FALSE)
  nwpgru

}
