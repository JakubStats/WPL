#--------------------------------------------------------------------------------------------------------
# The R-file given below contains functions/analysis used in the manuscript by Hwang, Heinze and Stoklosa 
# submitted to Biometrical Journal.
# 
# First, we give the necessary R-packages and functions to reproduce all figures, tables and simulations 
# presented in the manuscript.
#
# The analysis given below replicates our simulation studies and two examples. 
#
# The R-file presented below is split into 3 seperate parts (efficiency figures, two simulation studies 
# and two examples) where each part corresponds to figures and tables produced in the manuscript. The 
# results should be almost exact as set.seed() was used.
# 
# Further details regarding the simulation setup and model structure etc. are given in the manuscript
#
# Authors: Wen-Han Hwang (wenhan.hwang@gmail.com) and Jakub Stoklosa (j.stoklosa@unsw.edu.au).
#
# Please report any problems/suggestions to the authors listed above.
#
# These programs are meant to be used for non-commercial purposes only.
#--------------------------------------------------------------------------------------------------------

# Required R-Packages for all components given below (please install/load these first).

library(ggplot2)
library(VGAM)
library(reshape)
library(glmnet)
library(caret)
library(MASS)
library(mvtnorm)
library(RMark)
library(sandwich)
library(Rcapture)
library(numDeriv)
library(rootSolve)
library(nleqslv)
library(grid)

#--------------------------------------------------------------------------------------------------------

# **IMPORTANT**

# Below we give the required functions to generate, fit, plot and tabulate the results. Make sure these 
# functions are all stored and loaded first in your R console before running the sims and examples.

# Inverse logit function.

H <- function(a){1/(1+exp(-a))}

# Some simulation functions.

mean.na <- function(x) mean(x, na.rm = T)

# A function that calculates population size estimates (and standard errors).

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

# Efficiency comparison function for the homogeneous case.

mle.avar <- function(lam)
{
  q <- 1-exp(-lam)
  
  lam*q^2/(q - lam*exp(-lam))
}

pl.avar0 <- function(lam)
{
  p <- exp(-lam)
  q <- 1-exp(-lam)
  a <- (lam - q)/(lam*q)
  
  lam/a
}

pl.avar2 <- function(lam,K)
{
  p <- exp(-lam)
  q <- 1-exp(-lam)
  a <- (lam - q)/(lam*q)
  
  y <- 1:K
  g.m <- y-1-outer(y/(y + 1), lam)
  py.m <- t(outer(lam,y,"^")*p/q)
  py.m <- py.m/gamma(y + 1)
  eg2 <- apply(g.m^2*py.m,2,sum)
  
  eg2/(a^2)
}

binpl.avar <- function(p,k){
  k1 <- k+1
  q <- 1-p
  pi <- 1-q^k
  ey1.inv <- (1 - q^k1 - k1*p*q^k)/(k1*p*pi)
  a <- k-k1*ey1.inv
  y <- 1:k
  c <- (k*y - 1)/(y + 1)
  g.m <- y-1-outer(c, p)
  py.m0 <- (outer(p, y,  "^")*outer(q,k - y, "^"))/pi
  py.m <- t(py.m0)*choose(k, y)
  eg2 <- apply(g.m^2*py.m, 2, sum)
  
  eg2/(a^2)
}

binmle.avar <- function(p,k){
  q <- 1-p
  pi <- 1-q^k
  
  p*q*pi^2/(pi-k*p*q^(k-1))/k
}

binpl0.avar <- function(p,k){
  q <- 1-p
  pi <- 1-q^k
  
  p^2*q*pi/(k*p - pi)
}

# Multiple plot function.

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL){
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols*ceiling(numPlots / cols)), ncol = cols, nrow = ceiling(numPlots / cols))
  }
  
  if (numPlots == 1){
    print(plots[[1]])
  }else{
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  
  for(i in 1:numPlots){
      matchidx <- as.data.frame(which(layout==i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
  }
  }
}

#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# PART 1: Efficiency plots.
#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

# Positive-binomial models

# This should match figure 1 in the manscript.

Ns <- c(3, 5, 7, 10, 15, 20, 50, 75, 100)
sp <- seq(0.05, 0.95, by = 0.05)

ref <- list()

for(i in 1:length(Ns)){
  k <- Ns[i]
  
  bin.plavar <- binpl.avar(sp, k)
  bin.pl0avar <- binpl0.avar(sp, k)
  bin.mlavar <- binmle.avar(sp, k)
  
  bin.ref <- cbind(bin.mlavar / bin.pl0avar, bin.mlavar / bin.plavar)
  
  ref <- c(ref, list(bin.ref))
}

my.title1 <- bquote("m"==.(Ns[1]))
ref1 <- ref[[1]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_1 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.6, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m"==.(Ns[2]))
ref1 <- ref[[2]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_2 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.6, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m"==.(Ns[3]))
ref1 <- ref[[3]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_3 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.6, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m"==.(Ns[4]))
ref1 <- ref[[4]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_4 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.75, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m"==.(Ns[5]))
ref1 <- ref[[5]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_5 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.75, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m"==.(Ns[6]))
ref1 <- ref[[6]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_6 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.75, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m"==.(Ns[7]))
ref1 <- ref[[7]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_7 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.8, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m"==.(Ns[8]))
ref1 <- ref[[8]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_8 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.8, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m"==.(Ns[9]))
ref1 <- ref[[9]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_9 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.8, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

par(mfrow = c(3, 3), las=1)

multiplot(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, 
          layout = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                          ncol = 3, nrow = 3, byrow = T))

#--------------------------------------------------------------------------------------------------------

# Positive-Poisson models.

# This should match Figure 2 in the manscript.

slam <- seq(0.1, 20, by = 0.1)
plavar0 <- pl.avar0(slam)
plavar <- pl.avar2(slam, 50)
mlavar <- mle.avar(slam)
ef <- cbind(mlavar, plavar0, plavar)
ref <- ef[, 1]/ef[, 2:3]

# Plots.

par(mfrow = c(1, 1), las = 1)

colnames(ref) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref))
ref1 <- cbind(ref1,rep(slam, 2))
colnames(ref1) <- c("efficiency", "model", "lambda")

p_1 <- ggplot(ref1, aes(x = lambda, y = efficiency, colour = model, group = model)) +
  geom_line(size = 1, aes(linetype = model)) + xlab(expression(lambda)) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.8, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

multiplot(p_1, layout = matrix(c(1), ncol = 1, nrow = 1, byrow = T))

#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# PART 2: Simulations.
#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

# Simulation study 1a: Regression models.

# This should match Figure 3 and top of Table 1 in the manscript.

# Generate data, fit models and plot results here.

set.seed(312351)

simn <- 500 # No. of set simulation reps.

est.names <- c("PL", "WPL", "MLE")

ns <- c(50, 100, 200)

m <- 3 # Number of trials (or size) for a binomial model.

esttheta.MSE2 <- c()
esttheta.MSE_print <- c()

for(jjj in 1:4){
  esttheta.bias <- c()
  esttheta.MSE <- c()
  esttheta.MSE_n <- c()
  
  if(jjj == 1) theta <- c(-1, 1)
  if(jjj == 2) theta <- c(1.5, -0.5)
  if(jjj == 3) theta <- c(-1, 1)
  if(jjj == 4) theta <- c(1.5, -0.5)
  if(jjj == 5) theta <- c(-1, 1)
  if(jjj == 6) theta <- c(1.5, -0.5)
  
  for(n in ns){ 
    print(n)
    i <- 1
    esttheta1 <- matrix(NA, simn, length(est.names))
    esttheta.se1 <- matrix(NA, simn, length(est.names))
    esttheta.bias1 <- matrix(NA, simn, length(est.names))
    esttheta.MSE1 <- matrix(NA, simn, length(est.names))
    
    while(i <= simn){
      if (jjj==1 || jjj==2) x<-rnorm(n)
      
      if (jjj==3 || jjj==4){
        x <- rchisq(n, 3)
        x <- (x - mean(x))/sd(x)
        x[which(abs(x) > 3)] <- 3
      }
      if (jjj==5 || jjj==6){
        x <- runif(n, -3, 3)
        x <- (x - mean(x))/sd(x)
      }
      
      p <- H(theta[1] + theta[2]*x)
      y <- rposbinom(n, m, p)
      
      t_i<-c()
      
      for(jj in 1:length(y)){
        if (y[jj]==m){t1 <- 1}
        if (y[jj]!=m){
          pty <- choose(m - (1:m), y[jj] - 1)/choose(m, y[jj])
          t1 <- sample(1:m, 1, prob = pty)
        }
        
        t_i <- c(t_i, t1)
      }
      
      m.tilde <- m-t_i
      
      m.tilde.star <- m-(m + 1)/(y + 1)
      y.tilde <- y-1
      
      y.pl <- y.tilde/m.tilde
      y.pl[is.na(y.pl)] <- 0
      
      y.wpl <- y.tilde/m.tilde.star
      
      # Naive.    
      
      a1 <- glm(y / m ~ x, weights = rep(m, n), family = binomial)
      a1.se <- sqrt(vcov(a1)[2, 2])
      
      # PL.    
      
      a2 <- glm(y.pl ~ x, weights = m.tilde, family = binomial)
      a2.se <- sqrt(sandwich(a2)[2, 2])
      
      # WPL.
      
      a3 <- glm(y.wpl ~ x, weights = m.tilde.star, family = binomial)
      a3.se <- sqrt(sandwich(a3)[2, 2])
      
      # VGAM.
      
      dat1 <- data.frame(cbind(y, x))
      a4 <- try(vglm(cbind(y, m - y) ~ x, posbinomial, data = dat1), silent = TRUE)
      if(class(a4)=="try-error"){
        VGAM.coef <- NA
        a4.se <- NA
      }
      
      if(class(a4)!="try-error"){
        VGAM.coef <- coef(a4)[2]
        a4.se <- sqrt(vcov(a4)[2, 2])
      }
      
      esttheta1[i, ] <- c(coef(a2)[2], coef(a3)[2], VGAM.coef)
      esttheta.se1[i, ] <- c(a2.se, a3.se, a4.se)
      i <- i+1
    }
    
    esttheta.bias1 <- 100*apply(abs(esttheta1 - theta[2]) / theta[2], 2, mean.na)
    esttheta.bias <- rbind(esttheta.bias, esttheta.bias1)
    
    esttheta.MSE1 <- apply((esttheta1 - theta[2])^2, 2, mean.na)
    esttheta.MSE <- rbind(esttheta.MSE, esttheta.MSE1 / esttheta.MSE1[3])
    
    esttheta.MSE_n <- c(esttheta.MSE_n, c(round(esttheta.MSE1[1], digits = 3),
                                       9, round(esttheta.MSE1[2], digits = 3),
                                       9, round(esttheta.MSE1[3], digits = 3), 9))
  }
  
  esttheta.MSE_print <- rbind(esttheta.MSE_print, c(9, esttheta.MSE_n))
  
  esttheta.MSE <- esttheta.MSE[, -3]
  
  # Plots.
  
  colnames(esttheta.MSE) <- c("PL", "WPL")
  
  esttheta.MSE <- stack(as.data.frame(esttheta.MSE))
  
  if (jjj==1) beta_names <- "Scenario I"
  if (jjj==2) beta_names <- "Scenario II"
  if (jjj==3) beta_names <- "Scenario III"
  if (jjj==4) beta_names <- "Scenario IV"
  if (jjj==5) beta_names <- "Scenario V"
  if (jjj==6) beta_names <- "Scenario VI"
  
  esttheta.MSE <- cbind(esttheta.MSE, rep(ns, 2), rep(beta_names, 6))
  colnames(esttheta.MSE) <- c("Rel.MSE", "model", "n", "beta")
  
  esttheta.MSE2 <- rbind(esttheta.MSE2, esttheta.MSE)
}

colnames(esttheta.MSE_print) <- c("9", "PL", "9", "WPL", "9", "CL", "9", "PL", 
                                  "9", "WPL", "9", "CL", "9", "PL", "9", "WPL",
                                  "9", "CL", "9")

rownames(esttheta.MSE_print) <- c("Scenario I", "Scenario II", "Scenario III", "Scenario IV")

print(round(esttheta.MSE_print, digits = 3))

label <- "Simulation Study 1a"

df2 <- data.frame(Rel.MSE = c(1, 1), model = c("PL", "WPL"))

ggplot(data = esttheta.MSE2, aes(y = Rel.MSE, x = model, fill = model)) +
  scale_fill_brewer(palette = "Set2") + geom_col() + ggtitle(label, subtitle = NULL) +
  ylim(low = -0.01, high = 1.75) + facet_grid(beta ~ n, scales = "free") +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
  geom_bar(stat = "identity", colour = "black", position = "dodge") +
  geom_crossbar(data = df2, aes(x = factor(model), y = Rel.MSE, ymin = Rel.MSE, ymax = Rel.MSE),
                color = "black", linetype = "dotted")

#--------------------------------------------------------------------------------------------------------

# Simulation study 1b: Positive-Poisson simulations for log-linear regression.

# This should match Figure 4 and bottom of Table 1 in the manscript.

set.seed(20170922)

simn <- 500

est.names <- c("PL", "WPL", "MLE")

ns <- c(50, 100, 200)

esttheta.MSE2 <- c()
esttheta.MSE_print <- c()

for(jjj in 1:4)
  {
  esttheta.bias <- c()
  esttheta.MSE <- c()
  esttheta.MSE_n <- c()
  
  if(jjj==1) theta <- c(-1, 0.5)
  if(jjj==2) theta <- c(0.5, -1)
  if(jjj==3) theta <- c(-1, 0.5)
  if(jjj==4) theta <- c(0.5, -1)
  
  for(n in ns){
    print(n)
    i <- 1
    
    esttheta1 <- matrix(NA, simn, length(est.names))
    esttheta.se1 <- matrix(NA, simn, length(est.names))
    esttheta.bias1 <- matrix(NA, simn, length(est.names))
    esttheta.MSE1 <- matrix(NA, simn, length(est.names))
    
    while(i<=simn)
      {
      if(jjj==1 || jjj==2){x <- rnorm(n)}
      
      if(jjj==3 || jjj==4){
        x <- rchisq(n, 3)
        x <- (x - mean(x))/sd(x) 
        x[which(abs(x)>3)] <- 3
      }
      
      lambda <- exp(theta[1] + theta[2]*x)
      y <- rpospois(n, lambda)
      
      t_i <- c()
      
      for(jj in 1:length(y)){
        t1 <- min(runif(y[jj]))
        t_i <- c(t_i, t1)
      }
      
      t.tilde <- 1-t_i
      t.tilde.star <- y/(y + 1)
      y.tilde <- y-1
      
      # Naive.    
      
      a1 <- glm(y ~ x, family = poisson)
      a1.se <- sqrt(vcov(a1)[2, 2])
      
      # PL.
      
      a2 <- glm(y.tilde ~ x,offset = log(t.tilde), family = poisson)
      a2.se <- sqrt(sandwich(a2)[2, 2])
      
      # WPL.
      
      a3 <- glm(y.tilde ~ x, offset = log(t.tilde.star), family = poisson)
      a3.se <- sqrt(sandwich(a3)[2, 2])
      
      # VGAM.
      
      dat1 <- data.frame(cbind(y, x))
      a4 <- try(vglm(y ~ x, pospoisson, data = dat1), silent = TRUE)
      if (class(a4)=="try-error"){
        VGAM.coef <- NA
        a4.se <- NA
      }
      if (class(a4)!="try-error"){
        VGAM.coef <- coef(a4)[2]
        a4.se <- sqrt(vcov(a4)[2, 2])
      }
      
      esttheta1[i, ] <- c(coef(a2)[2], coef(a3)[2], VGAM.coef)
      esttheta.se1[i, ] <- c(a2.se, a3.se, a4.se)
      i <- i+1
      }
    
    esttheta.bias1 <- 100*apply(abs(esttheta1 - theta[2]) / theta[2], 2, mean.na)
    esttheta.bias <- rbind(esttheta.bias, esttheta.bias1)
    
    esttheta.MSE1 <- apply((esttheta1 - theta[2])^2, 2, mean.na)
    esttheta.MSE <- rbind(esttheta.MSE, esttheta.MSE1 / esttheta.MSE1[3])
    
    esttheta.MSE_n <- c(esttheta.MSE_n, c(round(esttheta.MSE1[1], digits = 3),
                                       9, round(esttheta.MSE1[2], digits = 3),
                                       9, round(esttheta.MSE1[3], digits = 3),
                                       9))
  }
  
  esttheta.MSE_print <- rbind(esttheta.MSE_print, c(9, esttheta.MSE_n))
  
  esttheta.MSE <- esttheta.MSE[, -3]
  
  # Plots.
  
  colnames(esttheta.MSE) <- c("PL", "WPL")
  
  esttheta.MSE <- stack(as.data.frame(esttheta.MSE))
  
  if (jjj==1) beta_names <- "Scenario I"
  if (jjj==2) beta_names <- "Scenario II"
  if (jjj==3) beta_names <- "Scenario III"
  if (jjj==4) beta_names <- "Scenario IV"
  if (jjj==5) beta_names <- "Scenario V"
  if (jjj==6) beta_names <- "Scenario VI"
  
  esttheta.MSE <- cbind(esttheta.MSE, rep(ns, 2), rep(beta_names, 6))
  colnames(esttheta.MSE) <- c("Rel.MSE", "model", "n", "beta")
  
  esttheta.MSE2 <- rbind(esttheta.MSE2, esttheta.MSE)
}

colnames(esttheta.MSE_print) <- c("9", "PL", "9", "WPL", "9", "CL", "9", "PL", "9", "WPL",
                                  "9", "CL", "9", "PL", "9", "WPL", "9", "CL", "9")
rownames(esttheta.MSE_print) <- c("Scenario I", "Scenario II", "Scenario III", "Scenario IV")

print(round(esttheta.MSE_print, digits = 3))

label <-"Simulation Study 1b"

df2 <- data.frame(Rel.MSE = c(1, 1), model = c("PL", "WPL"))

ggplot(data = esttheta.MSE2, aes(y = Rel.MSE, x = model, fill = model)) + scale_fill_brewer(palette = "Set2") + geom_col() +
  ggtitle(label,subtitle = NULL) + ylim(low = -0.01, high = 1.4) + facet_grid(beta ~ n, scales = "free") + 
  theme(legend.position = "none") + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
  geom_bar(stat = "identity", colour = "black", position = "dodge") + 
  geom_crossbar(data = df2, aes(x = factor(model), y = Rel.MSE, ymin = Rel.MSE, ymax = Rel.MSE), color = "black", linetype = "dotted")

#--------------------------------------------------------------------------------------------------------

# Simulation study 2: Variale selection and prediction.

# This should match Figures 5 and 6 of the manuscript (toggle with the family setting to get the desired figure).

set.seed(1)

simn <- 500

est.names <- c("MLE-Saturated", "WPL-AIC", "WPL-Elastic Net")

# Set the family here:

family <- "binomial"; m <- 10 # Figure 5.
#family <- "poisson"; m<-0  # Figure 6.

n <- 100 # Set the sample size.

if (family=="binomial") p <- 10; theta <- c(1, 3, 1.5, 0, 0, 2, 0, 0, 0, 0, 1)/3
if (family=="poisson") p <- 10; theta <- c(1, 3, 1.5, 0, 0, 2, 0, 0, 0, 0, 1)/3

start <- Sys.time()

SPR1 <- matrix(NA, simn, length(est.names))
DS1 <- matrix(NA, simn, length(est.names))
MSE1 <- matrix(NA, simn, length(est.names))
MAE1 <- matrix(NA, simn, length(est.names))

SPR2 <- matrix(NA, simn, length(est.names))
DS2 <- matrix(NA, simn, length(est.names))
MSE2 <- matrix(NA, simn, length(est.names))
MAE2 <- matrix(NA, simn, length(est.names))

i <- 1

while (i<=simn){
  X0 <- matrix(rep(1, n), ncol = 1)
  sigAR <- diag(p)
  sigAR <- 0.5^abs(row(sigAR) - col(sigAR))
  Xcov <- cbind(rmvnorm(n, mean = rep(0, p), sigAR))
  X <- as.data.frame(cbind(X0, Xcov))
  
  for(j in 1:p){
    colnames(X)[j+1] <- paste("x", j, sep = "")
  }
  colnames(X)[1] <- c("(Intercept)")
  X <- as.matrix(X)
  
  if (family=="binomial") pr <- H(drop(X%*%theta))
  if (family=="poisson")lambda <- exp(drop(X%*%theta))
  
  if (family=="binomial") y <- rposbinom(n, m, pr)
  if (family=="poisson") y <- rpospois(n, lambda)
  
  t_i <- c()
  
  for(jj in 1:length(y)){
    if (family=="binomial"){
      if (y[jj]==m) t1 <- 1
      if (y[jj]!=m){
        pty <- choose(m - (1:m), y[jj] - 1)/choose(m, y[jj])
        t1 <- sample(1:m, 1, prob = pty)
      }
    }
    
    if (family=="poisson") t1 <- min(runif(y[jj]))
    
    t_i <- c(t_i, t1)
  }
  
  if (family=="binomial"){
    m.tilde <- m - t_i
    m.tilde.star <- m-(m + 1)/(y + 1)
    y.tilde <- y-1
    y.wpl <- y.tilde/m.tilde.star
  }
  
  if (family=="poisson"){
    t.tilde.star <- y/(y + 1)
    y.tilde <- y-1
  }
  
  # Construct training/test data.
  
  train <- 1:80
  test <- 81:100
  n.tr <- length(train)
  
  n.te <- n-n.tr
  
  y.tr <- y[train]
  y.te <- y[test]
  
  X.tr <- X[train, ]
  X.te <- as.data.frame(X[test, ])
  
  if (family=="binomial"){
    y.wpl.tr <- y.wpl[train]
    y.wpl.te <- y.wpl[test]
    
    m.tilde.star.tr <- m.tilde.star[train]
    m.tilde.star.te <- m.tilde.star[test]
  }
  
  if (family=="poisson"){
    y.tilde.tr <- y.tilde[train]
    y.tilde.te <- y.tilde[test]
    
    t.tilde.star.tr <- t.tilde.star[train]
    t.tilde.star.te <- t.tilde.star[test]
  }  
  
  # Naive models.  
  
  # Zero-truncated models using vglm.
  
  dat1 <- data.frame(cbind(y.tr, X.tr))
  
  if (family=="binomial") a1 <- try(vglm(cbind(y.tr, m - y.tr) ~ X.tr - 1, posbinomial, data = dat1), silent = TRUE)
  if (family=="poisson") a1 <- try(vglm(y.tr ~ X.tr - 1, pospoisson, data = dat1), silent = TRUE)
  
  if(class(a1)=="try-error"){
    VGAM.SPR.tr <- NA
    VGAM.DS.tr <- NA
    VGAM.MSE.tr <- NA
    VGAM.MAE.tr <- NA
    
    VGAM.SPR.te <- NA
    VGAM.DS.te <- NA
    VGAM.MSE.te <- NA
    VGAM.MAE.te <- NA
  }
  
  if (class(a1)!="try-error"){
    dat2 <- data.frame(X.te)
    names(dat2) <- names(coef(a1))
    
    if (family=="binomial"){
      preds.tr0 <- c(H(coef(a1)%*%t(as.matrix(dat1[, -1]))))
      mu.tr <- m*preds.tr0
      pi.tr <- 1-(1 - mu.tr / m)^m
      y.hat.tr <- mu.tr/pi.tr
      var.hat.tr <- m*preds.tr0*(1 - preds.tr0)/pi.tr-m^2*preds.tr0^2*(1 - preds.tr0)^m/pi.tr^2
       
      VGAM.SPR.tr <- sum((y.tr - y.hat.tr)^2 / (var.hat.tr))/n.tr
      VGAM.DS.tr <- sum((y.tr - y.hat.tr)^2 / (var.hat.tr) + log(var.hat.tr))/n.tr
      VGAM.MSE.tr <- sum((y.tr - y.hat.tr)^2)/n.tr
      VGAM.MAE.tr <- sum(abs(y.tr - y.hat.tr))/n.tr
      
      preds.te0 <- c(H(coef(a1)%*%t(as.matrix(dat2))))
      mu.te <- m*preds.te0
      pi.te <- 1-(1 - mu.te / m)^m
      y.hat.te <- mu.te/pi.te
      var.hat.te <- m*preds.te0*(1 - preds.te0)/pi.te - m^2*preds.te0^2*(1 - preds.te0)^m/pi.te^2
      
      VGAM.SPR.te <- sum((y.te - y.hat.te)^2 / (var.hat.te))/n.te
      VGAM.DS.te <- sum((y.te - y.hat.te)^2 / (var.hat.te) + log(var.hat.te))/n.te
      VGAM.MSE.te <- sum((y.te - y.hat.te)^2)/n.te
      VGAM.MAE.te <- sum(abs(y.te - y.hat.te))/n.te
    }
    
    if (family=="poisson"){
      preds.tr0 <- c(exp(coef(a1)%*%t(as.matrix(dat1[,-1]))))
      preds.tr <- preds.tr0/(1 - exp(-preds.tr0))
      var.y.tr <- preds.tr-preds.tr^2*exp(-preds.tr0)
      
      VGAM.SPR.tr <- sum((y.tr - preds.tr)^2/var.y.tr)/n.tr
      VGAM.DS.tr <- sum((y.tr - preds.tr)^2/var.y.tr + log(var.y.tr))/n.tr
      VGAM.MSE.tr <- sum((y.tr - preds.tr)^2)/n.tr
      VGAM.MAE.tr <- sum(abs(y.tr - preds.tr))/n.tr
      
      preds.te0 <- c(exp(coef(a1)%*%t(as.matrix(dat2))))
      preds.te <- preds.te0/(1 - exp(-preds.te0))
      var.y.te <- preds.te-preds.te^2*exp(-preds.te0)
      
      VGAM.SPR.te <- sum((y.te - preds.te)^2 / var.y.te)/n.te
      VGAM.DS.te <- sum((y.te - preds.te)^2 / var.y.te + log(var.y.te))/n.te
      VGAM.MSE.te <- sum((y.te - preds.te)^2)/n.te
      VGAM.MAE.te <- sum(abs(y.te - preds.te))/n.te
    }
  }
  
  # WPL using StepAIC.
  
  if (family=="binomial") est <- try(glm(y.wpl.tr ~ X.tr - 1, weights = m.tilde.star.tr, family = binomial), silent = TRUE)
  if (family=="poisson"){
    offs <- log(t.tilde.star.tr)
    est <- try(glm(y.tilde.tr ~ X.tr - 1 + offset(offs), family = poisson), silent = TRUE)
  }
  
  if (class(est)[1]=="try-error"){
    WPL.AIC.SPR.tr <- NA
    WPL.AIC.DS.tr <- NA
    WPL.AIC.MSE.tr <- NA
    WPL.AIC.MAE.tr <- NA
    
    WPL.AIC.SPR.te <- NA
    WPL.AIC.DS.te <- NA
    WPL.AIC.MSE.te <- NA
    WPL.AIC.MAE.te <- NA
  }
  
  if (class(est)[1]!="try-error"){
    if (family=="binomial"){
      AIC.glm <- stepAIC(est, trace = FALSE, direction = c("backward"))$coefficients
      var.sel <- names(AIC.glm)
      dat2a <- X.tr[, which(names(coef(est))%in%var.sel)]
      preds.tr0 <- c(H(AIC.glm%*%t(as.matrix(dat2a))))
      mu.tr <- m*preds.tr0
      pi.tr <- 1-(1 - mu.tr / m)^m
      y.hat.tr <- mu.tr/pi.tr
      var.hat.tr<-m*preds.tr0*(1 - preds.tr0)/pi.tr-m^2*preds.tr0^2*(1 - preds.tr0)^m/pi.tr^2
      
      dat2b <- X.te[,which(names(coef(est))%in%var.sel)]
      preds.te0 <- c(H(AIC.glm%*%t(as.matrix(dat2b))))
      mu.te <- m*preds.te0
      pi.te <- 1-(1-mu.te / m)^m
      y.hat.te <- mu.te/pi.te
      var.hat.te <- m*preds.te0*(1 - preds.te0)/pi.te - m^2*preds.te0^2*(1 - preds.te0)^m/pi.te^2
      
      WPL.AIC.SPR.tr <- sum((y.tr - y.hat.tr)^2 / var.hat.tr)/n.tr
      WPL.AIC.DS.tr <- sum((y.tr - y.hat.tr)^2 / var.hat.tr + log(var.hat.tr))/n.tr
      WPL.AIC.MSE.tr <- sum((y.tr - y.hat.tr)^2)/n.tr
      WPL.AIC.MAE.tr <- sum((y.tr - y.hat.tr)^2)/n.tr
      
      WPL.AIC.SPR.te <- sum(((y.te - y.hat.te))^2 / (var.hat.te))/n.te
      WPL.AIC.DS.te <- sum(((y.te - y.hat.te))^2 / (var.hat.te) + log(var.hat.te))/n.te
      WPL.AIC.MSE.te <- sum((y.te - y.hat.te)^2)/n.te
      WPL.AIC.MAE.te <- sum((y.te - y.hat.te)^2)/n.te
    }
    
    if (family=="poisson"){
      AIC.glm <- stepAIC(est,trace = FALSE, direction = c("backward"))$coefficients
      var.sel <- names(AIC.glm)
      dat2a <- X.tr[, which(names(coef(est))%in%var.sel)]
      preds.tr0 <- c(exp(AIC.glm%*%t(as.matrix(dat2a))))
      preds.tr<-preds.tr0/(1 - exp(-preds.tr0))
      var.y.tr<-preds.tr-preds.tr^2*exp(-preds.tr0)
      
      dat2b <- X.te[, which(names(coef(est))%in%var.sel)]
      preds.te0 <- c(exp(AIC.glm%*%t(as.matrix(dat2b))))
      preds.te <- preds.te0/(1 - exp(-preds.te0))
      var.y.te <- preds.te-preds.te^2*exp(-preds.te0)
      
      WPL.AIC.SPR.tr <- sum((y.tr - preds.tr)^2 / var.y.tr)/n.tr
      WPL.AIC.DS.tr <- sum((y.tr - preds.tr)^2 / var.y.tr + log(var.y.tr))/n.tr
      WPL.AIC.MSE.tr <- sum((y.tr - preds.tr)^2)/n.tr
      WPL.AIC.MAE.tr <- sum(abs(y.tr - preds.tr))/n.tr
      
      WPL.AIC.SPR.te <- sum((y.te - preds.te)^2 / var.y.te)/n.te
      WPL.AIC.DS.te <- sum((y.te - preds.te)^2/var.y.te + log(var.y.te))/n.te
      WPL.AIC.MSE.te <- sum((y.te - preds.te)^2)/n.te
      WPL.AIC.MAE.te <- sum(abs(y.te - preds.te))/n.te
    }
  }
  
  # Elastic net variable selection.
  
  # WPL.
  
  if (family=="binomial"){
    Y.tr <- cbind(1 - y.wpl.tr, y.wpl.tr)
    est<-try(glmnet(as.matrix(X.tr[, -1]), Y.tr,
                    family = "binomial", weights = m.tilde.star.tr), silent = TRUE)
  }
  
  if (family=="poisson") est <- try(glmnet(X.tr[, -1], y.tilde.tr, family = "poisson", offset = log(t.tilde.star.tr)), silent = TRUE)
  if (class(est)[1]=="try-error"){
    WPL.GLMNET.SPR.tr <- NA
    WPL.GLMNET.DS.tr <- NA
    WPL.GLMNET.MSE.tr <- NA
    WPL.GLMNET.MAE.tr <- NA
    
    WPL.GLMNET.SPR.te <- NA
    WPL.GLMNET.MSE.te <- NA
    WPL.GLMNET.MAE.te <- NA
  }
  if (class(est)[1]!="try-error"){
    if (family=="binomial"){
      weight0.tr <- rep(1, n.tr)
      preds.tr0 <- predict(est,as.matrix(X.tr[, -1]),
                         s = cv.glmnet(as.matrix(X.tr[, -1]), Y.tr, family = "binomial", weights = m.tilde.star.tr)$lambda.min,
                         type = "response", newweights = weight0.tr)
      mu.tr <- m*preds.tr0
      pi.tr <- 1-(1 - mu.tr / m)^m
      y.hat.tr <- mu.tr/pi.tr
      var.hat.tr <- m*preds.tr0*(1 - preds.tr0)/pi.tr - m^2*preds.tr0^2*(1 - preds.tr0)^m/pi.tr^2
      
      weight0.te <- rep(1, n.te)
      preds.te0 <- predict(est,as.matrix(X.te[, -1]),
                         s = cv.glmnet(as.matrix(X.tr[, -1]), Y.tr, family = "binomial", weights = m.tilde.star.tr)$lambda.min,
                         type = "response", newweights = weight0.te)
      mu.te <- m*preds.te0
      pi.te <- 1-(1 - mu.te / m)^m
      y.hat.te <- mu.te/pi.te
      var.hat.te <- m*preds.te0*(1 - preds.te0)/pi.te-m^2*preds.te0^2*(1 - preds.te0)^m/pi.te^2
      
      WPL.GLMNET.SPR.tr <- sum((y.tr - y.hat.tr)^2 / var.hat.tr)/n.tr
      WPL.GLMNET.DS.tr <- sum((y.tr - y.hat.tr)^2 / var.hat.tr + log(var.hat.tr))/n.tr
      WPL.GLMNET.MSE.tr <- sum((y.tr - y.hat.tr)^2)/n.tr
      WPL.GLMNET.MAE.tr <- sum((y.tr - y.hat.tr)^2)/n.tr
      
      WPL.GLMNET.SPR.te <- sum(((y.te - y.hat.te)^2 / var.hat.te))/n.te
      WPL.GLMNET.DS.te <- sum((y.te - y.hat.te)^2 / var.hat.te + log(var.hat.te))/n.te
      WPL.GLMNET.MSE.te <- sum((y.te - y.hat.te)^2)/n.te
      WPL.GLMNET.MAE.te <- sum(abs(y.te - y.hat.te))/n.te
    }
    
    if (family=="poisson"){
      offset0.tr <- rep(0,n.tr)
      offset0.te <- rep(0,n.te)
      preds.tr0 <- predict(est,X.tr[, -1],
                         s = cv.glmnet(X.tr[, -1], y.tilde.tr, family = "poisson", offset = log(t.tilde.star.tr))$lambda.min,
                         type = "response", newoffset = offset0.tr)
      preds.te0 <- predict(est, as.matrix(X.te[, -1]),
                         s = cv.glmnet(X.tr[, -1], y.tilde.tr, family = "poisson", offset = log(t.tilde.star.tr))$lambda.min,
                         type = "response", newoffset = offset0.te)
      
      preds.tr <- preds.tr0/(1 - exp(-preds.tr0))
      preds.te <- preds.te0/(1 - exp(-preds.te0))
      
      var.y.tr <- preds.tr-preds.tr^2*exp(-preds.tr0)
      var.y.te <- preds.te-preds.te^2*exp(-preds.te0)
      
      WPL.GLMNET.SPR.tr <- sum((y.tr - preds.tr)^2 / var.y.tr)/n.tr
      WPL.GLMNET.DS.tr <- sum((y.tr - preds.tr)^2 / var.y.tr + log(var.y.tr))/n.tr
      WPL.GLMNET.MSE.tr <- sum((y.tr - preds.tr)^2)/n.tr
      WPL.GLMNET.MAE.tr <- sum(abs(y.tr - preds.tr))/n.tr
      
      WPL.GLMNET.SPR.te <- sum((y.te - preds.te)^2 / var.y.te)/n.te
      WPL.GLMNET.DS.te <- sum((y.te - preds.te)^2 / var.y.te + log(var.y.te))/n.te
      WPL.GLMNET.MSE.te <- sum((y.te - preds.te)^2)/n.te
      WPL.GLMNET.MAE.te <- sum(abs(y.te - preds.te))/n.te
    }
  }
  
  SPR1[i, ] <- c(VGAM.SPR.tr, WPL.AIC.SPR.tr, WPL.GLMNET.SPR.tr)
  DS1[i, ] <- c(VGAM.DS.tr, WPL.AIC.DS.tr, WPL.GLMNET.DS.tr)
  MSE1[i, ] <- c(VGAM.MSE.tr, WPL.AIC.MSE.tr, WPL.GLMNET.MSE.tr)
  MAE1[i, ] <- c(VGAM.MAE.tr, WPL.AIC.MAE.tr, WPL.GLMNET.MAE.tr)
  
  SPR2[i, ] <- c(VGAM.SPR.te, WPL.AIC.SPR.te, WPL.GLMNET.SPR.te)
  DS2[i, ] <- c(VGAM.DS.te, WPL.AIC.DS.te, WPL.GLMNET.DS.te)
  MSE2[i, ] <- c(VGAM.MSE.te, WPL.AIC.MSE.te, WPL.GLMNET.MSE.te)
  MAE2[i, ] <- c(VGAM.MAE.te, WPL.AIC.MAE.te, WPL.GLMNET.MAE.te)
  
  i <- i+1
}

DS.tr <- (DS1)
DS.te <- (DS2)
SPR.tr <- log(SPR1)
SPR.te <- log(SPR2)

colnames(DS.tr) = colnames(DS.te) = est.names
colnames(SPR.tr) = colnames(SPR.te) = est.names

meth.len <- length(est.names)

b1_mat <- stack(as.data.frame(DS.tr))
b2_mat <- stack(as.data.frame(SPR.tr))
b4_mat <- stack(as.data.frame(DS.te))
b5_mat <- stack(as.data.frame(SPR.te))

colnames(b1_mat) <- c("DS1", "method")
colnames(b2_mat) <- c("SPR1", "method")
colnames(b4_mat) <- c("DS2", "method")
colnames(b5_mat) <- c("SPR2", "method")

if (family=="poisson"){
  my.title1 <- c("Positive-Poisson (training) data")
  my.title2 <- c("Positive-Poisson (training) data")
  my.title4 <- c("Positive-Poisson (test) data")
  my.title5 <- c("Positive-Poisson (test) data")
}

if (family=="binomial"){
  my.title1 <- c("Positive-binomial (training) data")
  my.title2 <- c("Positive-binomial (training) data")
  my.title4 <- c("Positive-binomial (test) data")
  my.title5 <- c("Positive-binomial (test) data")
}

p1_gg <- ggplot(b1_mat, aes(method,DS1)) + labs(y = bquote("Dawid-Sebastiani score")) + 
  theme(axis.title.y = element_text(angle = 90, size = 16)) + theme(text = element_text(size = 14)) +
  geom_boxplot(aes(fill = method), outlier.colour = NA) + stat_boxplot(geom = 'errorbar') +
  ggtitle(my.title1) + scale_shape_identity() + theme(aspect.ratio = 1) +
  theme(legend.position = "none")

p2_gg <- ggplot(b2_mat, aes(method,SPR1)) + labs(y = bquote("Scaled Pearson residuals")) + 
  theme(axis.title.y = element_text(angle = 90, size = 16)) + theme(text = element_text(size = 14)) + 
  geom_boxplot(aes(fill = method), outlier.colour = NA) + stat_boxplot(geom = 'errorbar') +
  ggtitle(my.title2) + scale_shape_identity() + theme(aspect.ratio = 1) +
  theme(legend.position = "none")

p4_gg <- ggplot(b4_mat, aes(method,DS2)) + labs(y = bquote("Dawid-Sebastiani score")) +
  theme(axis.title.y = element_text(angle = 90, size = 16)) + theme(text = element_text(size = 14)) + 
  geom_boxplot(aes(fill = method),outlier.colour = NA) + stat_boxplot(geom = 'errorbar') + 
  ggtitle(my.title4) + scale_shape_identity() + theme(aspect.ratio = 1) + 
  theme(legend.position = "none")

p5_gg <- ggplot(b5_mat, aes(method,SPR2)) + labs(y = bquote("Scaled Pearson residuals")) + 
  theme(axis.title.y = element_text(angle = 90, size = 16)) + theme(text = element_text(size = 14)) + 
  geom_boxplot(aes(fill = method),outlier.colour = NA) + stat_boxplot(geom = 'errorbar') + 
  ggtitle(my.title5) + scale_shape_identity() + theme(aspect.ratio = 1) +
  theme(legend.position = "none")

multiplot(p2_gg, p5_gg, p1_gg, p4_gg, layout = matrix(c(1, 2, 3, 4), ncol = 2, nrow = 2, byrow = T))

#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# PART 3: Examples.
#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

# Example 1: Mt Little Higginbotham mountain pygmy possum data

# These should be the same as the first few rows of Table 3.

# These data can be obtained from Web Table 1.

tau <- 3 # No. of capture occasions.
D <- 62 # No. of uniquely caught possums.

# Observed frequencies:

y <- c(2, 3, 1, 3, 2, 2, 1, 3, 2, 2, 2, 1, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 1, 3, 
     1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 3, 2, 2, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 2, 1, 1, 2, 1, 3, 3, 1)

# First capture times (t_i) for each possum:

t1 <- c(1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 2, 2, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 2, 1, 3, 1, 1, 1, 2, 2, 2, 2, 2, 1, 3, 3, 3, 3, 3, 3, 
        2, 3, 3, 3, 3, 1, 1, 3, 1, 1, 2, 1, 1 ,1)

# Gender covariate:

x.obs <- c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0,
           1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 
           1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# Data summaries:

sum(y) # Total no. of captures.

mean(y) # Avergae capture rate.

sum(x.obs) # No. of males.

D - sum(x.obs) # No. of females.

fs <- table(y) 

f1 <- fs[1]
f2 <- fs[2]
f3 <- fs[3]

c(f1, f2, f3) # Frequency of individuals caught exactly x times.

# Use full data and fit linear logistic regression ($M_0/M_h$) models.

# Partial likelihood approach:=.

# Construct PL weights to feed into glm().

R <- y-1
h <- tau-t1
y.p <- R/h
y.p[is.na(y.p)] <- 0

est.PL_0 <- glm(y.p ~ 1, weights = h, family = binomial)
est.PL_const <- VarNhat.glm(est.PL_0, tau, y = y)

est.PL_1 <- glm(y.p ~ x.obs, weights = h, family = binomial)
est.PL_Mh <- VarNhat.glm(est.PL_1, tau, y = y)   

# Weighed partial likelihood approach.

# Construct WPL weights to feed into glm.

m.tilde.star <- tau-(tau + 1)/(y + 1)
h.wpl <- m.tilde.star
y.wpl <- R/h.wpl
y.wpl[is.na(y.wpl)] <- 0

est.WPL_0 <- glm(y.wpl ~ 1, weights = h.wpl, family = binomial)
est.WPL_const <- VarNhat.glm(est.WPL_0, tau, y = y)

est.WPL_1 <- glm(y.wpl ~ x.obs, weights = h.wpl, family = binomial)
est.WPL_Mh <- VarNhat.glm(est.WPL_1, tau, y = y) 

# Combine results and display them.

N_ests <- matrix(round(as.numeric(rbind(c(est.PL_const$Nhat, est.PL_const$Se.Nhat),
                                      c(est.WPL_const$Nhat, est.WPL_const$Se.Nhat),
                                      c(est.PL_Mh$Nhat, est.PL_Mh$Se.Nhat),
                                      c(est.WPL_Mh$Nhat, est.WPL_Mh$Se.Nhat))),
                                      digits = 2), ncol = 2)

rownames(N_ests) <- c("PL", "WPL", "PL-h", "WPL-h")
colnames(N_ests) <- c("N_hat", "S.E.(N_hat)")
round(N_ests, digits = 2)

#--------------------------------------------------------------------------------------------------------

# Example 2: Variable selection using GLMNET.

# Should be similar the second column of Table 4.

# The 1987/88 US National Medical Expenditure Survey (NMES) count data were obtained from: 
# https://www.jstatsoft.org/article/view/v027i08
 
# Load data and extract all variables.

load(file = "DebTrivedi.rda") # Load the data.

dt <- DebTrivedi[,c(1, 5, 6, 8, 9, 11:19)]
dt[, 5] <- as.numeric(dt[, 5])-1
dt[, 7] <- as.numeric(dt[, 7])-1
dt[, 8] <- as.numeric(dt[, 8])-1
dt[, 9] <- as.numeric(dt[, 9])-1
dt[, 12] <- as.numeric(dt[, 12])-1
dt[, 13] <- as.numeric(dt[, 13])-1
dt[, 14] <- as.numeric(dt[, 14])-1

# Remove all zero counts from data to create artificial zero-truncated data.

dt2 <- dt[-which(dt$ofp==0), ]

y <- dt2$ofp
n <- length(y)

X <- cbind(rep(1, n), dt2[, -1])
colnames(X)[1] <- "(Intercept)"

# Fit models and apply model selection (AIC and GLMNET).

# Construct WPL weights to feed into glm() and glmnet().

t.tilde.star <- y/(y + 1)
y.tilde <- y-1

# Fit models here.

dat2 <- data.frame(cbind(y.tilde, X))

offs <- log(t.tilde.star)

mod2 <- glm(y.tilde ~ emer + hosp + numchron + adldiff + age + black + 
              gender + married + school + faminc + employed + privins + medicaid,
            offset = offs, family = poisson, data = dat2)

AIC.glm2 <- stepAIC(mod2, trace = FALSE)$coefficients

mod3 <- glmnet(as.matrix(X), y.tilde, family = "poisson", offset = log(t.tilde.star))

s <- cv.glmnet(as.matrix(X), y.tilde, family = "poisson", offset = log(t.tilde.star))$lambda.min

tmp_coeffs <- coef(mod3, s = s)

mod3.coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

# Combine results and display them.

AIC.glm2 

# These will be slightly different from the last column of Table 4 because glmnet() uses 
# cross-validation to select lambda, thus the data is randomly split and will consist of 
# different training/test sets for each fit.

t(mod3.coef)

