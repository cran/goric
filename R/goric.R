goric <-
function(object, ..., iter=100000){
  if (!inherits(object, "orlm")) stop("object needs to be of class orlm")
  if (iter < 1) stop("No of iterations < 1")
  objlist <- list(object, ...)
  isorlm <- sapply(objlist, function(x) inherits(x, "orlm"))
  orlmlist <- objlist[isorlm]  
  Call <- match.call()
  Call$iter <- NULL
  names(orlmlist) <- as.character(Call[-1L])[isorlm]
  loglik <- sapply(orlmlist, function(x) x$logLik)
  penalty <- sapply(orlmlist, function(x) goric_penalty(x, iter=iter))
  goric <- 2*(loglik - penalty)
  goric_weights <- exp(goric/2) / sum(exp(goric/2))
  data.frame(loglik, penalty, goric=-1*goric, goric_weights=round(goric_weights,3))
}

