# soft-thresholding operator
S <- function(z, lambda) 
{
  (z - lambda) * (z > lambda) + 
    (z + lambda) * (z < -lambda) + 
    0 * (abs(z) <= lambda)
}


# SCAD update function for CD
scad.update <- function(z, lambda, a = 3.7)
{
  if (abs(z) < 2 * lambda){
    S(z, lambda)
  } else if ((2 * lambda <= abs(z)) & (abs(z) <= a * lambda)){
    S(z, a * lambda / (a - 1)) / (1 - 1 / (a - 1))
  } else z
}



# coordinate decent algorithm
ncv <- function(x, y, lambda, type = "lasso", a = 3.7, init = rep(0, p), max.iter = 100, eps = 1.0e-8)
{
  n <- length(y)
  p <- ncol(x)
  
  # marginal standardization of x
  x <- scale(x)
  m <- attr(x, "scaled:center")
  s <- attr(x, "scaled:scale")
  
  # centering of y
  my <- mean(y)
  y <- (y - my)
  
  # initialize beta
  beta <- init
  
  # residual
  r <- (y - x %*% beta)
  
  if (type == "lasso") 
  {
    update.ft <- function(x, lambda) lasso.update(x, lambda)
  } else if (type == "scad") {
    update.ft <- function(x, lambda) scad.update(x, lambda, a)
  } else if (type == "mcp") {
    update.ft <- function(x, lambda) mcp.update(x, lambda, a)
  } else stop("type should be lasso, scad, or mcp!")
  
  # start update
  for (t in 1:max.iter)
  {
    new.beta <- beta
    for (j in 1:p)
    {
      xj <- 1/n * crossprod(x[,j],  r) + beta[j]
      new.beta[j] <- update.ft(xj, lambda/s[j]) 
      r <- r - (new.beta[j] - beta[j]) * x[,j]  
    }
    if (max(abs(beta - new.beta)) < eps) break
    beta <- new.beta
  }
  
  # transform back
  beta <- beta / s[j]
  beta0 <- my - m %*% beta
  
  index <- which(abs(beta) > eps)
  beta.info <- beta[index]
  
  obj = list(intercept = beta0,
             beta = beta.info,
             index = index)
}

#LQA algorithm


scad_lqa <- function(x, y, lambda, a = 3.7, init = rep(1, p), max.iter = 100, eps = 1.0e-8)
{
  n <- length(y)
  p <- ncol(x)
  
  # marginal standardization of x
  x <- scale(x)
  m <- attr(x, "scaled:center")
  s <- attr(x, "scaled:scale")
  
  # centering of y
  my <- mean(y)
  y <- (y - my)
  yhat <- y
  
  # initialize beta
  beta <- init
  
  # penalty
  
  dp <- function(x){ifelse(x<=lambda,lambda,lambda*(a*lambda-x)/((a-1)*lambda))}
 
  
  # start update
  for (t in 1:max.iter)
  {
    new.beta <- beta
    for (j in 1:p)
    {
      if (new.beta[j] == 0) next
      wj <- dp(abs(new.beta[j]))/(2*abs(new.beta[j]))
      xj <- ((t(x)%*%y)/n)[j]
      new.beta[j] <- xj/(1+2*wj)
      if(abs(new.beta[j]) < eps){
        new.beta[j] = 0
        next;}
    }
    if (max(abs(beta - new.beta)) < eps) break
    beta <- new.beta
  }
  
  # transform back
  
  beta <- beta / s[j]
  beta0 <- my - m %*% beta
  
  index <- which(abs(beta) > 1.0e-2)
  beta.info <- beta[index]
  
  obj = list(intercept = beta0,
             beta = beta.info,
             index = index)
}



# LLA algorithm
scad_lla <- function(x, y, lambda, a = 3.7, init = rep(0, p), max.iter = 100, eps = 1.0e-8)
{
  n <- length(y)
  p <- ncol(x)
  
  # marginal standardization of x
  x <- scale(x)
  m <- attr(x, "scaled:center")
  s <- attr(x, "scaled:scale")
  
  # centering of y
  my <- mean(y)
  y <- (y - my)
  
  # initialize beta
  beta <- init
  
  # penalty
  dp <- function(x){ifelse(x<=lambda,lambda,lambda*(a*lambda-x)/((a-1)*lambda))}
  update.ft <- function(x, lambda) scad.update(x, lambda, a)
  
  # start update
  for (t in 1:max.iter)
  {
    new.beta <- beta
    for (j in 1:p)
    {
      wj <- dp(abs(new.beta[j]))
      xj <- ((t(x)%*%y)/n)[j]
      new.beta[j] <- update.ft(xj, wj/s[j]) 
    }
    if (max(abs(beta - new.beta)) < eps) break
    beta <- new.beta
  }
  
  # transform back
  beta <- beta / s[j]
  beta0 <- my - m %*% beta
  
  index <- which(abs(beta) > 1.0e-2)
  beta.info <- beta[index]
  
  obj = list(intercept = beta0,
             beta = beta.info,
             index = index)
}


# One step algorithm
scad_one <- function(x, y, lambda, a = 3.7, eps = 1.0e-8)
{
  n <- length(y)
  p <- ncol(x)
  
  # marginal standardization of x
  x <- scale(x)
  m <- attr(x, "scaled:center")
  s <- attr(x, "scaled:scale")
  
  # centering of y
  my <- mean(y)
  y <- (y - my)
  
  # initialize beta to MLE
  beta <- apply(x,2,mean)
  
  # penalty
  dp <- function(x){ifelse(x<=lambda,lambda,lambda*(a*lambda-x)/((a-1)*lambda))}
  update.ft <- function(x, lambda) scad.update(x, lambda, a)
  
  # start update
     new.beta <- beta
    for (j in 1:p)
    {
      wj <- dp(abs(new.beta[j]))
      xj <- ((t(x)%*%y)/n)[j]
      new.beta[j] <- update.ft(xj, wj/s[j]) 
    beta <- new.beta
  }
  
  # transform back
  beta <- beta / s[j]
  beta0 <- my - m %*% beta
  
  index <- which(abs(beta) > 1.0e-2)
  beta.info <- beta[index]
  
  obj = list(intercept = beta0,
             beta = beta.info,
             index = index)
}




set.seed(1)
n <- 100
p <- 5

d <- 2
beta <- c(rep(1, d), rep(0, p-d))

x <- matrix(rnorm(n*p), n, p)
e <- rnorm(n)
y <- x %*% beta + e

lambda <- 0.2


# scad cd
obj1 <- ncv(x, y, lambda, type = "scad")
# scad lqa
obj2 <- scad_lqa(x, y, lambda, a = 3.7)
# scad lla
obj3 <- scad_lla(x, y, lambda, a = 3.7)
# scad one-step
obj4 <- scad_one(x, y, lambda, a = 3.7)


print(obj1)
print(obj2)
print(obj3)
print(obj4)

