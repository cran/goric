goric_penalty <-
function(object, iter=100000, mc.cores=1){
  require(quadprog)
  require(mvtnorm)
  if (!inherits(object, "orlm")) stop("object needs to be of class orlm")
  if (all(object$constr == 0) & object$nec == 0){
    penalty <- length(object$coefficients) + 1
  } else {
    if (iter < 1) stop("No of iterations < 1")
    Sigma <- object$sigma
    x <- object$X
    invW <- kronecker(solve(Sigma), t(x) %*% x)
    W <- solve(invW)
    Z <- rmvnorm(n=iter, mean=rep(0, ncol(W)), sigma=W)
    Dmat <- 2*invW
    Amat <- t(object$constr)
    bvec <- object$rhs
    nec <- object$nec
    if (mc.cores == 1){
      nact <- apply(Z, 1, function(z){
        dvec <- 2*(z %*% invW)
        QP <- solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=nec)
        length(QP$iact)
      })
    } else {
      require(parallel)
      cl <- makeCluster(mc.cores)
      clusterExport(cl, c("solve.QP"))
      nact <- parRapply(cl, Z, function(z, invW,Dmat,Amat,bvec,nec){
        dvec <- 2*(z %*% invW)
        QP <- solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=nec)
        length(QP$iact)
      }, invW=invW, Dmat=Dmat, Amat=Amat, bvec=bvec, nec=nec)
      stopCluster(cl)
    }
    dimsol <- ncol(W) - nact
    LP <- sapply(1:(ncol(W)+1), function(x) sum(x == (dimsol+1)))/iter
    penalty <- 1 + sum((1:ncol(W))*LP[-1])
  }
  return(penalty)
}

