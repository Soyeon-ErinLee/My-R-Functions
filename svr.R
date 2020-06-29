library(kernlab) 
# Support Vector Regression

svr <- function(x, y, e = 1, lambda = 1) {
  # x: predictor
  # y: response
  # e: epsilon of the e-loss
  # lambda: regularization parameter
  
  n = length(y)
  H <- x %*% t(x)
  c <- rep(e,n)-y
  b <- 0
  r <- 0
  A <- matrix(1,1,n)
  l <- rep(-1/lambda,n)
  u <- rep((1/lambda),n)
  obj <- ipop(c,H,A,b,l,u,r)
  theta <- obj@primal
  beta <- c(crossprod(theta,x))
  temp1<-temp2<-matrix(0)
  eps <- 1.0e-7
  sv.index <- which(abs(theta)>=eps)
  ytemp <- y[sv.index]
  x <- as.matrix(x)
  xtemp <- as.matrix(x[sv.index,])
  ntemp <- length(sv.index)
  if(ntemp!=0){
    for(i in 1:ntemp){
      ifelse(ytemp[i]<0, temp1[i] <- ytemp[i]-xtemp[i,]%*%beta+e, temp2[i]<-ytemp[i]-xtemp[i,]%*%beta-e)
      beta0 <- (sum(temp1)+sum(temp2))/ntemp
    }
  } else{
    beta0 <- 0}
  est <- c(beta0, beta)
  # beta0 for intercept, and beta for the coefficient vector
  return(est)
}
