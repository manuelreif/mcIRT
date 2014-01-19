deviance.nelm <-
function(object, ...)
{
  nme  <- length(RESnlm$erg_distr$mean_est) - 1
  nva  <- RESnlm$ctrl$sigmaest *(length(RESnlm$erg_distr$sig_est) -1)
  npar <- ncol(RESnlm$reshOBJ$Qmat) + nme + nva - length(RESnlm$ctrl$Clist)
  
structure(2*object$last_mstep$value, df=npar)
}
