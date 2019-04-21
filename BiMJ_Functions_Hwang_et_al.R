#--------------------------------------------------------------------------------------------------------
# BiMJ_Functions_Hwang_et_al.R
#--------------------------------------------------------------------------------------------------------

# Required R-Packages for all components given below (please install/load these first).

library(ggplot2)
library(VGAM)
library(reshape)
library(glmnet)
library(caret)
library(MASS)
library(mvtnorm)
library(sandwich)
library(Rcapture)
library(numDeriv)
library(rootSolve)
library(nleqslv)
library(grid)

# Required functions to generate, fit, plot and tabulate the results. 

# Inverse logit function.

H <- function(a){1/(1+exp(-a))}

# Some simulation functions.

mean.na <- function(x) mean(x, na.rm = T)

# A function that calculates population size estimates (and standard errors) given a glm model.

VarNhat.glm <- function(m, tau, y = NULL){
  
  # The "m" denotes a glm model object, tau and y are as defined as above.
  
  X <- model.matrix(m)
  beta <- coef(m)
  P <- c(1 / (1 + exp(-X%*%beta)))
  Pi <- 1 - (1 - P)^tau
  Nhat <- sum(Pi^(-1)) # Population size (Horvitz-Thompson) estimator.
  
  # Standard error estimator using the "sandwich" package. Below we give the
  # standard error estimator for Nhat, see Huggins (1989) for further details.
  
  var.beta <- sandwich(m) # Variance estimates for model regression coefficients (beta).
  gdash.beta <- t(X)%*%(Pi^(-2)*(1 - P)^tau*tau*P)
  varA<-sum((1 - Pi)/Pi^2)
  varB <- (t(gdash.beta)%*%var.beta)%*%gdash.beta
  varN <- as.vector(varA + varB)
  Se.Nhat <- sqrt(varN)
  
  return(list(Se.beta = sqrt(diag(var.beta)), Nhat = Nhat, Se.Nhat = Se.Nhat))
}

# Conditional likelihood function for the positive-binomial model.
# Populations size (and stadnard errors) are also estimated. 

par.mle <- function(betain, X, tau, y, SE = T){
  options(warn = -1)
  
  n <- nrow(X)
  
  est.fn <- function(beta){
    pr <- H(X%*%beta)
    pi <- 1-(1-pr)^tau
    L <- y*log(pr/(1-pr))+tau*log(1-pr)-log(pi)
    
    -sum(L)
  }
  
  est.out <- optim(betain, est.fn, control = list(maxit = 10000))
  beta <- est.out$par
  
  pr <- c(H(X%*%beta))
  pi <- 1-(1-pr)^tau
  N.hat <- sum(1/pi)
  
  if(SE==F){
    sd.beta <- NULL
    sd.Nhat <- NULL
  }
  
  # Variance calculation.
  
  if(SE==T){
    vy <- tau*pr*(1-pr)/pi-(tau*pr)^2*(1-pi)/pi/pi
    dG <- t(X)%*%(X*vy)
    idG <- solve(dG)
    sd.beta <- sqrt(diag(idG))
    
    dpi<-tau*pr*(1-pi)/pi/pi
    dNhat<-apply(X*dpi, 2, sum)
    sd.Nhat<-sqrt(sum((1-pi)/pi/pi)+t(dNhat)%*%idG%*%dNhat)
  }
  
  return(list(beta = beta, N.hat = N.hat, sd.Nhat = sd.Nhat, sd.b = sd.beta, pr = pr))
}

# Efficiency comparison function for the homogeneous case.

mle.avar <- function(lam){
  q <- 1-exp(-lam)
  
  lam*q^2/(q - lam*exp(-lam))
}

pl.avar0 <- function(lam){
  p <- exp(-lam)
  q <- 1-exp(-lam)
  a <- (lam - q)/(lam*q)
  
  lam/a
}

pl.avar2 <- function(lam, K){
  p <- exp(-lam)
  q <- 1-exp(-lam)
  a <- (lam - q)/(lam*q)
  
  y <- 1:K
  g.m <- y-1-outer(y/(y + 1), lam)
  py.m <- t(outer(lam,y,"^")*p/q)
  py.m <- py.m/gamma(y + 1)
  eg2 <- apply(g.m^2*py.m, 2, sum)
  
  eg2/(a^2)
}

binpl.avar <- function(p, k){
  k1 <- k + 1 
  q <- 1- p
  pi <- 1 - q^k
  ey1.inv <- (1 - q^k1 - k1*p*q^k)/(k1*p*pi)
  a <- k - k1*ey1.inv
  y <- 1:k
  c <- (k*y - 1)/(y + 1)
  g.m <- y-1-outer(c, p)
  py.m0 <- (outer(p, y,  "^")*outer(q, k - y, "^"))/pi
  py.m <- t(py.m0)*choose(k, y)
  eg2 <- apply(g.m^2*py.m, 2, sum)
  
  eg2/(a^2)
}

binmle.avar <- function(p,k){
  q <- 1-p
  pi <- 1-q^k
  
  p*q*pi^2/(pi-k*p*q^(k-1))/k
}

binpl0.avar <- function(p, k){
  q <- 1-p
  pi <- 1-q^k
  
  p^2*q*pi/(k*p - pi)
}

# Multiple plot function.

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL){
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  
  if (is.null(layout)){
    layout <- matrix(seq(1, cols*ceiling(numPlots / cols)), ncol = cols, nrow = ceiling(numPlots / cols))
  }
  
  if (numPlots  ==  1){
    print(plots[[1]])
  }else{
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for(i in 1:numPlots){
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}