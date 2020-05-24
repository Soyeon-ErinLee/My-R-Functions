# Kernel Quadratic Regression

library(kernlab) 
my_kqr <- function(x,y,n,C = 1,tau){
    H <- x %*% t(x)
    c <- -y
    b <- 0
    r <- 0
    A <- matrix(1,1,n)
    l <- rep((tau-1),n)
    u <- rep(tau,n)
    obj <- ipop(c,H,A,b,l,u,r)
    theta <-  obj@primal
    beta <- c(crossprod(theta,x))
    sv.index<-which((tau-1)*C < theta & theta < tau*C)
    temp <- y[sv.index]-x[sv.index,]%*%beta
    beta0 <- mean(temp)
    yhat <- beta0 + x%*%beta
    resid <- y-yhat
    loss <- tau*resid - resid*(resid<0)
    return(list(beta0=beta0, beta=beta,resid=resid,loss=loss))
  }

set.seed(1)
n <- 100
C <- 1
tau <- 0.5
p <- 5
beta <- rep(1,p)
x <- matrix(rnorm(n*p), n, p)
e<-rnorm(n,0,0)
y <- x %*% beta + e

my_kqr(x,y,n,C,tau)

par(mfrow=c(1,1))
plot(my_kqr(x,y,n,C,tau=0.20)$resid,my_kqr(x,y,n,C,tau=0.20)$loss,xlab = "residual", ylab = "loss",xlim=c(-6e-11,6e-11),ylim=c(0,3e-11),pch=1)
par(new=TRUE)
plot(my_kqr(x,y,n,C,tau=0.50)$resid,my_kqr(x,y,n,C,tau=0.50)$loss,xlab = "residual", ylab = "loss",xlim=c(-6e-11,6e-11),ylim=c(0,3e-11),pch=2)
par(new=TRUE)
plot(my_kqr(x,y,n,C,tau=0.75)$resid,my_kqr(x,y,n,C,tau=0.75)$loss,xlab = "residual", ylab = "loss",xlim=c(-6e-11,6e-11),ylim=c(0,3e-11),pch=3)
legend('top',legend = paste("tau = ", c(0.20,0.50,0.75), sep=''), pch=c(1,2,3))
