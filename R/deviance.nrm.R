deviance.nrm <-
function(object, ...)
{
  nme  <- length(object$erg_distr$mean_est) - 1
  nva  <- length(object$erg_distr$sig_est) -1
  npar <- ncol(object$reshOBJ$Qmat) + nme + nva  
structure(2*object$last_mstep$value, df=npar)
}
