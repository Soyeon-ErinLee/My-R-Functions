# Soft-Tresholding function
S <- function(z, lambda) {
  (z - lambda) * (z > lambda) + 
    (z + lambda) * (z < -lambda) + 
    0 * (abs(z) <= lambda)
}

# Coordinate decent algorithm for LASSO
cd.lasso <- function(x, y, lambda) 
{
  z <- scale(x)
  m <- attr(z, "scaled:center")
  s <- attr(z, "scaled:scale")
  u <- (y - mean(y))
  beta <- coef(lm(u ~ z - 1))
  r <- u - z %*% beta
  
  for (iter in 1:100) {
    new.beta <- beta
    for (j in 1:p) {
      temp <- beta[j] + crossprod(z[,j], r)/n
      new.beta[j] <- S(temp, lambda/s[j]) 
      r <- r - (new.beta[j] - beta[j]) * z[,j]
    }
    delta <- max(abs(new.beta - beta))
    if (delta < 1.0e-3) break
    beta <- new.beta
  }
  beta <- new.beta/s 
  intercept <- mean(y) - crossprod(beta, m)
  
  obj<-list(intercept=intercept, beta=beta)
  return(obj)
}


# Pathwise Coordinate decent algorithm for LASSO
pathcd.lasso <- function(x, y, x.test, y.test, grid = 100)
{
  n <- length(y)
  p <- ncol(x)
  u <- (y - mean(y))
  z <- scale(x)

  lambda <- seq(max(abs(crossprod(z, y))), 0, length = grid)
  beta <- matrix(0, grid, p)
  beta0 <- rep(0, grid)
  mse <- rep(0,grid)
  for (i in 1:grid){
    cd <- cd.lasso(x, y, lambda[i])
    beta0[i] <- cd$intercept
    beta[i,] <- cd$beta
    yhat <- beta0[i] + x.test %*% beta[i,]
    mse[i] <- sum((y.test - yhat)^2)
  }
  index <- which.min(mse)
  obj <- list(intercept = beta0[index], beta=beta[index,], lambda=lambda[index])  
  return(obj)
}  

train <- matrix(scan("train.txt"), 500, 51)
test <- matrix(scan("test.txt"), 500, 51)
x <- train[,-51]
y <- train[,51]
x.test <- test[,-51]
y.test <- test[,51]

pathcd.lasso(x,y,x.test,y.test)
