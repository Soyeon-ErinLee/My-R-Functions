# Monte Carlo Markov Chain Cauchy Using Metropolis Hastings Algorithm

set.seed(1)
n.sim <- 100000
X <- c(rt(1,1)) # t-distribution becomes Cauchy distribution when df=1
ac <- NULL # acceptance rate
ac[1] <- 1

# Markov Chain
for (i in 2:n.sim) {
  Y <- rnorm(1)
  rho <- dt(Y,1)*dnorm(X[i-1]) / (dt(X[i-1],1)*dnorm(Y))
  ac[i] <- runif(1) < rho
  X[i] <- X[i-1] + (Y - X[i-1]) * ac[i]
}

cat("acceptance rate = ", mean(ac), "\b")


hist(X, probability = T, breaks = 50, main = "MH algoritghm",xlim = c(-10,10))
lines(seq(-10,10,length=1000), dt(seq(-10,10,length=1000),1), col = 2, lty = 2)

obj <- ks.test(jitter(X), rt(n.sim, 1))
print(obj)


# The chain based on the normal proposal is consistently off the true value
