goric <-
function(object, ..., iter=100000, mc.cores=1){
  if (!inherits(object, "orlm") & !inherits(object, "list")) stop("object needs to be of class orlm or a list of orlm objects")
  if (iter < 1) stop("No of iterations < 1")
  if (inherits(object, "orlm")) objlist <- list(object, ...) else objlist <- object
  isorlm <- sapply(objlist, function(x) inherits(x, "orlm"))
  orlmlist <- objlist[isorlm]  
  Call <- match.call()
  Call$iter <- NULL
  Call$mc.cores <- NULL
  if (inherits(object, "orlm")) names(orlmlist) <- as.character(Call[-1L])[isorlm]
  loglik <- sapply(orlmlist, function(x) x$logLik)
  penalty <- sapply(orlmlist, function(x) goric_penalty(x, iter=iter, mc.cores=mc.cores))
  goric <- -2*(loglik - penalty)
  delta <- goric - min(goric)
  goric_weights <- exp(-delta/2) / sum(exp(-delta/2))
  data.frame(loglik, penalty, goric=goric, delta=round(delta,3), goric_weights=round(goric_weights,3))
}

