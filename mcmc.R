# Monte Carlo Markov Chain

mcmc <- function(x, y, n.samples = 10000,  init = rep(0, p+1), step = rep(0.3,p+1)) {
  # x: n * p predictor matrix
  # y: response vector in the Poisson regression
  # n.samples: number of posterior samples to be obtained.
  # init, step: other factors required in MCMC
  
  p <- ncol(x)
  x1 <- cbind(rep(1, n), x)
  
  post.beta <- matrix(0, n.samples, p+1)
  prior.m <- 0 
  prior.s <- 1000
  post.beta[1,] <- beta <- init
  eta <- x1%*%beta 
  mu <- exp(eta)
  
  log.like <- sum(-mu+y*log(mu))
  
  for (iter in 1:n.samples){
    beta.new <- beta 
    for (j in 1:(p+1))
    {
      beta.new[j] <- beta[j] + rnorm(1, 0, step[j]) 
      eta.new <- x1%*%beta.new
      mu.new <- exp(eta.new)
      
      log.prior <- dnorm(beta[j], prior.m, prior.s, log = T)
      log.prior.new <- dnorm(beta.new[j], prior.m, prior.s, log = T)
      log.like.new <- sum(-mu.new+y*log(mu.new))
      
      temp <- exp((log.like.new + log.prior.new) - (log.like + log.prior))
      rho <- min(1, temp)
      
      if (runif(1) < rho) {
        beta[j] <- beta.new[j]
        log.like  <- log.like.new
        eta <- x1 %*% beta
        mu <- exp(eta)
      }
    }
    post.beta[iter,] <- beta
  }
  samples <- post.beta
  return(samples) 
}
