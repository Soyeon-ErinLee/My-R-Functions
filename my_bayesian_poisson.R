
mcmc.poisson <- function(x, y, init, n.sample = 10000, step = 0.3){
  post.beta <- matrix(0, n.sample, p)
  ac.ratio <- rep(0, n.sample)
  
  prior.m <- 0
  prior.s <- 1000 # for vague prior
  
  # intialize
  post.beta[1,] <- beta <- init
  eta <- x %*% beta
  mu <- exp(eta)
  
  log.prior <- sum(dnorm(beta, prior.m, prior.s, log = T))
  log.like  <- sum(-mu + y * log(mu))
  
  iter <- 2
  for (iter in 1:n.sample){
    
    # candidate
    beta.new <- beta + rnorm(p, 0, step) 
    eta.new <- x %*% beta.new
    mu.new  <- exp(eta.new)
    
    # prior
    log.prior.new <- sum(dnorm(beta.new, prior.m, prior.s, log = T))
    
    # liklihood
    log.like.new <- sum(-mu.new + y * log(mu.new))
    
    # ratio
    temp <- exp((log.like.new + log.prior.new) - (log.like + log.prior))
    rho <- min(1, temp)
    
    if (runif(1) < rho) {
      ac.ratio[iter] <- 1
      beta <- beta.new
      log.prior <- log.prior.new
      log.like  <- log.like.new
      eta <- x %*% beta
      mu <- exp(eta)
    }
    post.beta[iter,] <- beta
  }
  
  obj <- list(posterior = post.beta, acpt.ratio = mean(ac.ratio))
  return(obj)
}


# Bayesian Logistic Regression
set.seed(1)

n <- 100 # sample size
p <- 1   # predictor dimension

x <- matrix(runif(n*p), n, p) # generate predictor

true.beta <- 0 # true beta
true.mu <- exp(true.beta * x)  
y <- rpois(n, true.mu)         # generate reponse


init <- rep(5, p)


# mle 
obj.mle <- glm(y ~ x - 1, family = "poisson")
mle <- coef(obj.mle)
cat("ML estimate =", mle, "\n")
print("CI based on MLE:")
print(confint(obj.mle))

# MH for Bayesian Poisson
obj <- mcmc.poisson(x, y, init = init, n.sample = 10000, step = 0.3)
posterior <- obj$posterior[-(1:2000),]
est <- mean(posterior)
cat("accptance ratio =",obj$acpt.ratio, "\n")
cat("Bayesian estimate =", est, "\n")

cr <- t(apply(as.matrix(posterior), 2, quantile, c(0.025, 0.975))) # Credible Region
print("CR for posteiors:")
print(cr)

#comparison
library(MCMCpack)
obj2 <- summary(MCMCpoisson(y~x-1))
obj2$statistics
obj2$quantiles


# trace plot
par(mfrow = c(1,1))
for (k in 1:p){
  plot(obj$posterior[,k], type = "l", col = "gray", ylab = "posterior", xlab = "", main = "trace plot")
  abline(h = est[k], col = "blue", lwd = 2)
  abline(h = mle[k], col = "green", lty = 3, lwd = 2)
  abline(h = true.beta[k], col = "red", lty = 2, lwd = 2)
  legend("topright", c("mle", "Bayes", "True"), col = c(3,4,2), lty = c(3, 1, 2))
}

