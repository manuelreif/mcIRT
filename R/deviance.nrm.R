deviance.nrm <-
function(object, ...)
{
  ### number of parameters
  nme  <- length(RESnrm$erg_distr$mean_est) - 1
  #RESnrm$ctrl$sigmaest
  nva  <- RESnrm$ctrl$sigmaest * (length(RESnrm$erg_distr$sig_est) -1)
  npar <- ncol(RESnrm$reshOBJ$Qmat) + nme + nva  - length(RESnrm$ctrl$Clist)
  
structure(2*object$last_mstep$value, df=npar)
}
