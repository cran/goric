goric <-
function(object, ..., iter=100000, mc.cores=1){
UseMethod("goric")
}


goric.orlm <-
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
  data.frame(loglik, penalty, goric=goric, goric_weights=round(goric_weights,3))
}


goric.orgls <-
function(object, ..., iter=100000, mc.cores=1){
  if (!inherits(object, "orgls") & !inherits(object, "list")) stop("object needs to be of class orgls or a list of orgls objects")
  if (iter < 1) stop("No of iterations < 1")
  if (inherits(object, "orgls")) objlist <- list(object, ...) else objlist <- object
  isorgls <- sapply(objlist, function(x) inherits(x, "orgls"))
  orlmlist <- objlist[isorgls]  
  Call <- match.call()
  Call$iter <- NULL
  Call$mc.cores <- NULL
  if (inherits(object, "orgls")) names(orlmlist) <- as.character(Call[-1L])[isorgls]
  loglik <- sapply(orlmlist, function(x) x$logLik)
  ep <- sapply(orlmlist, function(x) x$extrap)
  penalty <- sapply(orlmlist, function(x) goric_penalty(x, iter=iter, mc.cores=mc.cores))
  goric <- -2*(loglik - penalty - ep)
  delta <- goric - min(goric)
  goric_weights <- exp(-delta/2) / sum(exp(-delta/2))
  data.frame(loglik, penalty, vcdf=ep, goric=goric, goric_weights=round(goric_weights,3))
}

goric.list <- function(object, ..., iter=100000, mc.cores=1){
  if (all(sapply(object, class) == "orlm")) out <- goric.orlm(object, iter=iter, mc.cores=mc.cores)
  if (all(sapply(object, class) == "orgls")) out <- goric.orgls(object, iter=iter, mc.cores=mc.cores)
  return(out)
}
