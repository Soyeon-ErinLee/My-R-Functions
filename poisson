
# Poisson beta estimation

my.poisson <- function(X, y, init = rep(1, ncol(X)), maxiter = 100, eps = 1.0e-8){
  beta <- init
  iter <- 0
  while (iter < maxiter) {
    mu <- exp(as.vector(X%*%beta))
    W <- diag(mu)
    Xtil <- sqrt(W)%*%X
    ytil <- sqrt(W)%*%(X%*%beta + diag(1/mu)%*%(y-mu))
    qr.obj <- qr(Xtil)
    new.beta <- backsolve(qr.obj$qr,qr.qty(qr.obj,ytil))
    if (sum((new.beta - beta)^2)/sum(beta^2) < eps) break
    beta <- new.beta
    iter <- iter + 1
  }
  if (iter == maxiter) warning("Algorithm may not be converged!")
  obj <- list(beta = new.beta, iteration = iter)
  return(obj)
}

 
