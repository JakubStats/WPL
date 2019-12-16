#--------------------------------------------------------------------------------------------------------
# BiMJ_Simulations_Hwang_et_al.R
#--------------------------------------------------------------------------------------------------------

rm(list = ls())

# Make sure all functions stored and loaded first in your R-console.

source("BiMJ_Functions_Hwang_et_al.R")

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

for(jjj in 1:4) {
  esttheta.bias <- c()
  esttheta.MSE <- c()
  esttheta.MSE_n <- c()
  
  if(jjj  ==  1) theta <- c(-1, 1)
  if(jjj  ==  2) theta <- c(1.5, -0.5)
  if(jjj  ==  3) theta <- c(-1, 1)
  if(jjj  ==  4) theta <- c(1.5, -0.5)
  if(jjj  ==  5) theta <- c(-1, 1)
  if(jjj  ==  6) theta <- c(1.5, -0.5)
  
  for(n in ns) { 
    print(n)
    i <- 1
    esttheta1 <- matrix(NA, simn, length(est.names))
    esttheta.se1 <- matrix(NA, simn, length(est.names))
    esttheta.bias1 <- matrix(NA, simn, length(est.names))
    esttheta.MSE1 <- matrix(NA, simn, length(est.names))
    
    while(i <= simn) {
      if (jjj == 1 || jjj == 2) x<-rnorm(n)
      
      if (jjj == 3 || jjj == 4) {
        x <- rchisq(n, 3)
        x <- (x - mean(x))/sd(x)
        x[which(abs(x) > 3)] <- 3
      }
      if (jjj == 5 || jjj == 6) {
        x <- runif(n, -3, 3)
        x <- (x - mean(x))/sd(x)
      }
      
      p <- H(theta[1] + theta[2]*x)
      y <- rposbinom(n, m, p)
      
      t_i<-c()
      
      for(jj in 1:length(y)) {
        if (y[jj] == m) {t1 <- 1}
        if (y[jj] != m) {
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
      if(class(a4) == "try-error") {
        VGAM.coef <- NA
        a4.se <- NA
      }
      
      if(class(a4) != "try-error") {
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
  
  if (jjj == 1) beta_names <- "Scenario I"
  if (jjj == 2) beta_names <- "Scenario II"
  if (jjj == 3) beta_names <- "Scenario III"
  if (jjj == 4) beta_names <- "Scenario IV"
  if (jjj == 5) beta_names <- "Scenario V"
  if (jjj == 6) beta_names <- "Scenario VI"
  
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

rm(list = ls())

# Make sure all functions stored and loaded first in your R-console.

source("BiMJ_Functions_Hwang_et_al.R")

# This should match Figure 4 and bottom of Table 1 in the manscript.

set.seed(20170922)

simn <- 500

est.names <- c("PL", "WPL", "MLE")

ns <- c(50, 100, 200)

esttheta.MSE2 <- c()
esttheta.MSE_print <- c()

for(jjj in 1:4) {
  esttheta.bias <- c()
  esttheta.MSE <- c()
  esttheta.MSE_n <- c()
  
  if(jjj == 1) theta <- c(-1, 0.5)
  if(jjj == 2) theta <- c(0.5, -1)
  if(jjj == 3) theta <- c(-1, 0.5)
  if(jjj == 4) theta <- c(0.5, -1)
  
  for(n in ns) {
    print(n)
    i <- 1
    
    esttheta1 <- matrix(NA, simn, length(est.names))
    esttheta.se1 <- matrix(NA, simn, length(est.names))
    esttheta.bias1 <- matrix(NA, simn, length(est.names))
    esttheta.MSE1 <- matrix(NA, simn, length(est.names))
    
    while(i<=simn) {
      if(jjj == 1 | jjj == 2) {x <- rnorm(n)}
      
      if(jjj == 3 | jjj == 4) {
        x <- rchisq(n, 3)
        x <- (x - mean(x))/sd(x) 
        x[which(abs(x)>3)] <- 3
      }
      
      lambda <- exp(theta[1] + theta[2]*x)
      y <- rpospois(n, lambda)
      
      t_i <- c()
      
      for(jj in 1:length(y)) {
        t1 <- min(runif(y[jj]))
        t_i <- c(t_i, t1)
      }
      
      t.tilde <- 1 - t_i
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
      if (class(a4) == "try-error") {
        VGAM.coef <- NA
        a4.se <- NA
      }
      if (class(a4) != "try-error") {
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
  
  if (jjj == 1) beta_names <- "Scenario I"
  if (jjj == 2) beta_names <- "Scenario II"
  if (jjj == 3) beta_names <- "Scenario III"
  if (jjj == 4) beta_names <- "Scenario IV"
  if (jjj == 5) beta_names <- "Scenario V"
  if (jjj == 6) beta_names <- "Scenario VI"
  
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

ggplot(data = esttheta.MSE2, aes(y = Rel.MSE, x = model, fill = model)) + 
  scale_fill_brewer(palette = "Set2") + geom_col() +
  ggtitle(label,subtitle = NULL) + ylim(low = -0.01, high = 1.4) + 
  facet_grid(beta ~ n, scales = "free") + 
  theme(legend.position = "none") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
  geom_bar(stat = "identity", colour = "black", position = "dodge") + 
  geom_crossbar(data = df2, aes(x = factor(model), y = Rel.MSE, ymin = Rel.MSE, ymax = Rel.MSE), 
                color = "black", linetype = "dotted")

#--------------------------------------------------------------------------------------------------------

# Simulation study 2: Variale selection and prediction.

rm(list = ls())

# Make sure all functions stored and loaded first in your R-console.

source("BiMJ_Functions_Hwang_et_al.R")

# These should match Figures 5 and 6 of the manuscript (make sure to toggle with 
# the family setting below to get the desired figure).

set.seed(1)

simn <- 500

est.names <- c("MLE-Saturated", "WPL-AIC", "WPL-Elastic Net")

# Set the family here:

#family <- "binomial"; m <- 10 # Figure 5.
family <- "poisson"; m<-0  # Figure 6.

n <- 100 # Set the sample size.

if (family == "binomial") p <- 10; theta <- c(1, 3, 1.5, 0, 0, 2, 0, 0, 0, 0, 1)/3
if (family == "poisson") p <- 10; theta <- c(1, 3, 1.5, 0, 0, 2, 0, 0, 0, 0, 1)/3

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

while (i<=simn) {
  X0 <- matrix(rep(1, n), ncol = 1)
  sigAR <- diag(p)
  sigAR <- 0.5^abs(row(sigAR) - col(sigAR))
  Xcov <- cbind(rmvnorm(n, mean = rep(0, p), sigAR))
  X <- as.data.frame(cbind(X0, Xcov))
  
  for(j in 1:p) {
    colnames(X)[j+1] <- paste("x", j, sep = "")
  }
  colnames(X)[1] <- c("(Intercept)")
  X <- as.matrix(X)
  
  if (family == "binomial") pr <- H(drop(X%*%theta))
  if (family == "poisson") lambda <- exp(drop(X%*%theta))
  
  if (family == "binomial") y <- rposbinom(n, m, pr)
  if (family == "poisson") y <- rpospois(n, lambda)
  
  t_i <- c()
  
  for(jj in 1:length(y)) {
    if (family == "binomial") {
      if (y[jj] == m) t1 <- 1
      if (y[jj] != m) {
        pty <- choose(m - (1:m), y[jj] - 1)/choose(m, y[jj])
        t1 <- sample(1:m, 1, prob = pty)
      }
    }
    
    if (family == "poisson") t1 <- min(runif(y[jj]))
    
    t_i <- c(t_i, t1)
  }
  
  if (family == "binomial") {
    m.tilde <- m - t_i
    m.tilde.star <- m-(m + 1)/(y + 1)
    y.tilde <- y-1
    y.wpl <- y.tilde/m.tilde.star
  }
  
  if (family == "poisson") {
    t.tilde.star <- y/(y + 1)
    y.tilde <- y - 1
  }
  
  # Construct training/test data.
  
  train <- 1:80
  test <- 81:100
  n.tr <- length(train)
  
  n.te <- n - n.tr
  
  y.tr <- y[train]
  y.te <- y[test]
  
  X.tr <- X[train, ]
  X.te <- as.data.frame(X[test, ])
  
  if (family == "binomial") {
    y.wpl.tr <- y.wpl[train]
    y.wpl.te <- y.wpl[test]
    
    m.tilde.star.tr <- m.tilde.star[train]
    m.tilde.star.te <- m.tilde.star[test]
  }
  
  if (family == "poisson") {
    y.tilde.tr <- y.tilde[train]
    y.tilde.te <- y.tilde[test]
    
    t.tilde.star.tr <- t.tilde.star[train]
    t.tilde.star.te <- t.tilde.star[test]
  }  
  
  # Naive models.  
  
  # Zero-truncated models using vglm.
  
  dat1 <- data.frame(cbind(y.tr, X.tr))
  
  if (family == "binomial") a1 <- try(vglm(cbind(y.tr, m - y.tr) ~ X.tr - 1, 
                                           posbinomial, data = dat1), silent = TRUE)
  if (family == "poisson") a1 <- try(vglm(y.tr ~ X.tr - 1, 
                                          pospoisson, data = dat1), silent = TRUE)
  
  if(class(a1) == "try-error") {
    VGAM.SPR.tr <- NA
    VGAM.DS.tr <- NA
    VGAM.MSE.tr <- NA
    VGAM.MAE.tr <- NA
    
    VGAM.SPR.te <- NA
    VGAM.DS.te <- NA
    VGAM.MSE.te <- NA
    VGAM.MAE.te <- NA
  }
  
  if (class(a1) != "try-error") {
    dat2 <- data.frame(X.te)
    names(dat2) <- names(coef(a1))
    
    if (family == "binomial") {
      preds.tr0 <- c(H(coef(a1)%*%t(as.matrix(dat1[, -1]))))
      mu.tr <- m*preds.tr0
      pi.tr <- 1 - (1 - mu.tr / m)^m
      y.hat.tr <- mu.tr/pi.tr
      var.hat.tr <- m*preds.tr0*(1 - preds.tr0)/pi.tr-m^2*preds.tr0^2*(1 - preds.tr0)^m/pi.tr^2
      
      VGAM.SPR.tr <- sum((y.tr - y.hat.tr)^2 / (var.hat.tr))/n.tr
      VGAM.DS.tr <- sum((y.tr - y.hat.tr)^2 / (var.hat.tr) + log(var.hat.tr))/n.tr
      VGAM.MSE.tr <- sum((y.tr - y.hat.tr)^2)/n.tr
      VGAM.MAE.tr <- sum(abs(y.tr - y.hat.tr))/n.tr
      
      preds.te0 <- c(H(coef(a1)%*%t(as.matrix(dat2))))
      mu.te <- m*preds.te0
      pi.te <- 1 - (1 - mu.te / m)^m
      y.hat.te <- mu.te/pi.te
      var.hat.te <- m*preds.te0*(1 - preds.te0)/pi.te - m^2*preds.te0^2*(1 - preds.te0)^m/pi.te^2
      
      VGAM.SPR.te <- sum((y.te - y.hat.te)^2 / (var.hat.te))/n.te
      VGAM.DS.te <- sum((y.te - y.hat.te)^2 / (var.hat.te) + log(var.hat.te))/n.te
      VGAM.MSE.te <- sum((y.te - y.hat.te)^2)/n.te
      VGAM.MAE.te <- sum(abs(y.te - y.hat.te))/n.te
    }
    
    if (family == "poisson") {
      preds.tr0 <- c(exp(coef(a1)%*%t(as.matrix(dat1[, -1]))))
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
  
  if (family == "binomial") est <- try(glm(y.wpl.tr ~ X.tr - 1, weights = m.tilde.star.tr, 
                                           family = binomial), silent = TRUE)
  if (family == "poisson") {
    offs <- log(t.tilde.star.tr)
    est <- try(glm(y.tilde.tr ~ X.tr - 1 + offset(offs), family = poisson), silent = TRUE)
  }
  
  if (class(est)[1] == "try-error") {
    WPL.AIC.SPR.tr <- NA
    WPL.AIC.DS.tr <- NA
    WPL.AIC.MSE.tr <- NA
    WPL.AIC.MAE.tr <- NA
    
    WPL.AIC.SPR.te <- NA
    WPL.AIC.DS.te <- NA
    WPL.AIC.MSE.te <- NA
    WPL.AIC.MAE.te <- NA
  }
  
  if (class(est)[1] != "try-error") {
    if (family == "binomial") {
      AIC.glm <- stepAIC(est, trace = FALSE, direction = c("backward"))$coefficients
      var.sel <- names(AIC.glm)
      dat2a <- X.tr[, which(names(coef(est))%in%var.sel)]
      preds.tr0 <- c(H(AIC.glm%*%t(as.matrix(dat2a))))
      mu.tr <- m*preds.tr0
      pi.tr <- 1 - (1 - mu.tr / m)^m
      y.hat.tr <- mu.tr/pi.tr
      var.hat.tr<-m*preds.tr0*(1 - preds.tr0)/pi.tr - m^2*preds.tr0^2*(1 - preds.tr0)^m/pi.tr^2
      
      dat2b <- X.te[,which(names(coef(est))%in%var.sel)]
      preds.te0 <- c(H(AIC.glm%*%t(as.matrix(dat2b))))
      mu.te <- m*preds.te0
      pi.te <- 1 - (1 - mu.te / m)^m
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
    
    if (family == "poisson") {
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
  
  if (family == "binomial") {
    Y.tr <- cbind(1 - y.wpl.tr, y.wpl.tr)
    est<-try(glmnet(as.matrix(X.tr[, -1]), Y.tr,
                    family = "binomial", weights = m.tilde.star.tr), silent = TRUE)
  }
  
  if (family == "poisson") est <- try(glmnet(X.tr[, -1], y.tilde.tr, family = "poisson",
                                             offset = log(t.tilde.star.tr)), silent = TRUE)
  if (class(est)[1] == "try-error") {
    WPL.GLMNET.SPR.tr <- NA
    WPL.GLMNET.DS.tr <- NA
    WPL.GLMNET.MSE.tr <- NA
    WPL.GLMNET.MAE.tr <- NA
    
    WPL.GLMNET.SPR.te <- NA
    WPL.GLMNET.MSE.te <- NA
    WPL.GLMNET.MAE.te <- NA
  }
  if (class(est)[1] != "try-error") {
    if (family == "binomial") {
      weight0.tr <- rep(1, n.tr)
      preds.tr0 <- predict(est,as.matrix(X.tr[, -1]),
                           s = cv.glmnet(as.matrix(X.tr[, -1]), Y.tr, family = "binomial", 
                                         weights = m.tilde.star.tr)$lambda.min,
                           type = "response", newweights = weight0.tr)
      mu.tr <- m*preds.tr0
      pi.tr <- 1 - (1 - mu.tr / m)^m
      y.hat.tr <- mu.tr/pi.tr
      var.hat.tr <- m*preds.tr0*(1 - preds.tr0)/pi.tr - m^2*preds.tr0^2*(1 - preds.tr0)^m/pi.tr^2
      
      weight0.te <- rep(1, n.te)
      preds.te0 <- predict(est,as.matrix(X.te[, -1]),
                           s = cv.glmnet(as.matrix(X.tr[, -1]), Y.tr, family = "binomial", 
                                         weights = m.tilde.star.tr)$lambda.min,
                           type = "response", newweights = weight0.te)
      mu.te <- m*preds.te0
      pi.te <- 1 - (1 - mu.te / m)^m
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
    
    if (family == "poisson") {
      offset0.tr <- rep(0,n.tr)
      offset0.te <- rep(0,n.te)
      preds.tr0 <- predict(est,X.tr[, -1],
                           s = cv.glmnet(X.tr[, -1], y.tilde.tr, family = "poisson", 
                                         offset = log(t.tilde.star.tr))$lambda.min,
                           type = "response", newoffset = offset0.tr)
      preds.te0 <- predict(est, as.matrix(X.te[, -1]),
                           s = cv.glmnet(X.tr[, -1], y.tilde.tr, family = "poisson", 
                                         offset = log(t.tilde.star.tr))$lambda.min,
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

if (family == "poisson") {
  my.title1 <- c("Positive-Poisson (training) data")
  my.title2 <- c("Positive-Poisson (training) data")
  my.title4 <- c("Positive-Poisson (test) data")
  my.title5 <- c("Positive-Poisson (test) data")
}

if (family == "binomial") {
  my.title1 <- c("Positive-binomial (training) data")
  my.title2 <- c("Positive-binomial (training) data")
  my.title4 <- c("Positive-binomial (test) data")
  my.title5 <- c("Positive-binomial (test) data")
}

p1_gg <- ggplot(b1_mat, aes(method,DS1)) + labs(y = bquote("Dawid-Sebastiani score")) + 
  theme(axis.title.y = element_text(angle = 90, size = 16)) + theme(text = element_text(size = 10)) +
  geom_boxplot(aes(fill = method), outlier.colour = NA) + stat_boxplot(geom = 'errorbar') +
  ggtitle(my.title1) + scale_shape_identity() + theme(aspect.ratio = 1) +
  theme(legend.position = "none")

p2_gg <- ggplot(b2_mat, aes(method,SPR1)) + labs(y = bquote("Scaled Pearson residuals")) + 
  theme(axis.title.y = element_text(angle = 90, size = 16)) + theme(text = element_text(size = 10)) + 
  geom_boxplot(aes(fill = method), outlier.colour = NA) + stat_boxplot(geom = 'errorbar') +
  ggtitle(my.title2) + scale_shape_identity() + theme(aspect.ratio = 1) +
  theme(legend.position = "none")

p4_gg <- ggplot(b4_mat, aes(method,DS2)) + labs(y = bquote("Dawid-Sebastiani score")) +
  theme(axis.title.y = element_text(angle = 90, size = 16)) + theme(text = element_text(size = 10)) + 
  geom_boxplot(aes(fill = method),outlier.colour = NA) + stat_boxplot(geom = 'errorbar') + 
  ggtitle(my.title4) + scale_shape_identity() + theme(aspect.ratio = 1) + 
  theme(legend.position = "none")

p5_gg <- ggplot(b5_mat, aes(method,SPR2)) + labs(y = bquote("Scaled Pearson residuals")) + 
  theme(axis.title.y = element_text(angle = 90, size = 16)) + theme(text = element_text(size = 10)) + 
  geom_boxplot(aes(fill = method),outlier.colour = NA) + stat_boxplot(geom = 'errorbar') + 
  ggtitle(my.title5) + scale_shape_identity() + theme(aspect.ratio = 1) +
  theme(legend.position = "none")

multiplot(p2_gg, p5_gg, p1_gg, p4_gg, layout = matrix(c(1, 2, 3, 4), ncol = 2, nrow = 2, byrow = T))
