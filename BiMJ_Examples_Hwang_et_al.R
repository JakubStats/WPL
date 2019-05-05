#--------------------------------------------------------------------------------------------------------
# BiMJ_Simulations_Hwang_et_al.R
#--------------------------------------------------------------------------------------------------------

rm(list = ls())

# Make sure all functions stored and loaded first in your R-console.

source("BiMJ_Functions_Hwang_et_al.R")

# Example 1: Mt Little Higginbotham mountain pygmy possum data.

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

cap.hist <- read.table("cap.hist.txt") # Load the capture history data.

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

# Closed form estimates for males/females separately. 

# "N" denotes the population size estimate throughout.

# These results should be the same as Table 2.

female <- which(x.obs==0)
male <- which(x.obs==1)

# Partial likelihood approach:

# Construct PL weights to feed into glm().

t0 <- t(apply(cap.hist, 1, cumsum))
t0a <- apply(t0==0, 1, sum)+1
t1 <- t0a[y>0]

R <- y-1
h <- tau-t1
y.p <- R/h
y.p[is.na(y.p)] <- 0

p_PL_female <- sum(R[female])/sum(h[female])
p_PL_male <- sum(R[male])/sum(h[male])

pi_PL_female <- rep((1-(1-p_PL_female)^tau), length(y[female]))
pi_PL_male <- rep((1-(1-p_PL_male)^tau), length(y[male]))

N_PL_female <- sum(1/pi_PL_female)
N_PL_male <- sum(1/pi_PL_male)

N_PL_tot <- N_PL_male + N_PL_female

est.PL0_female <- glm(y.p[female] ~ 1, weights = h[female], family = binomial)
est.PL_female <- VarNhat.glm(est.PL0_female, tau, y = y[female])    
est.PL0_male <- glm(y.p[male] ~ 1, weights = h[male], family = binomial)
est.PL_male <- VarNhat.glm(est.PL0_male, tau, y = y[male])

N_PL_tot.SE <- sqrt(est.PL_female$Se.Nhat^2 + est.PL_male$Se.Nhat^2)

# Weighed partial likelihood approach:

# Construct WPL weights to feed into glm.

m.tilde.star <- tau - (tau+1)/(y+1)
h.wpl <- m.tilde.star
y.wpl <- R/h.wpl
y.wpl[is.na(y.wpl)] <- 0

p_WPL_female <- sum(R[female])/sum(m.tilde.star[female])
p_WPL_male <- sum(R[male])/sum(m.tilde.star[male])

pi_WPL_female <- rep((1-(1-p_WPL_female)^tau), length(y[female]))
pi_WPL_male <- rep((1-(1-p_WPL_male)^tau), length(y[male]))

N_WPL_female <- sum(1/pi_WPL_female)
N_WPL_male <- sum(1/pi_WPL_male)

N_WPL_tot <- N_WPL_male + N_WPL_female

est.WPL0_female <- glm(y.wpl[female] ~ 1, weights = h.wpl[female], family = binomial)
est.WPL_female <- VarNhat.glm(est.WPL0_female, tau, y = y[female])    
est.WPL0_male <- glm(y.wpl[male] ~ 1, weights = h.wpl[male], family = binomial)
est.WPL_male <- VarNhat.glm(est.WPL0_male, tau, y = y[male])

N_WPL_tot.SE <- sqrt(est.WPL_female$Se.Nhat^2 + est.WPL_male$Se.Nhat^2)

# Conditional likelihood (CL).

Z_female <- matrix(c(rep(1, length(x.obs[female]))), ncol = 1) 
est.CL_female <- par.mle(betain = 1, Z_female, tau, y = y[female]) 

Z_male <- matrix(c(rep(1, length(x.obs[male]))), ncol = 1) 
est.CL_male <- par.mle(betain = 1, Z_male, tau, y = y[male]) 

N_CL_female <- est.CL_female$N.hat
N_CL_male <- est.CL_male$N.hat

N_CL_tot <- N_CL_male + N_CL_female

N_CL_tot.SE <- sqrt(est.CL_female$sd.Nhat^2 + est.CL_male$sd.Nhat^2)

N_ests <- rbind(c(p_PL_female, sqrt(sandwich(est.PL0_female))[1]*est.PL0_female$fitted[1]*(1-est.PL0_female$fitted[1]), 
                  p_PL_male, sqrt(sandwich(est.PL0_male))[1]*est.PL0_male$fitted[1]*(1-est.PL0_male$fitted[1]), NA, NA), 
            c(p_WPL_female, sqrt(sandwich(est.WPL0_female))[1]*est.WPL0_female$fitted[1]*(1-est.WPL0_female$fitted[1]), 
              p_WPL_male, sqrt(sandwich(est.WPL0_male))[1]*est.WPL0_male$fitted[1]*(1-est.WPL0_male$fitted[1]), NA, NA), 
            c(H(est.CL_female$beta), (est.CL_female$sd.b)*(est.CL_female$pr[1]*(1-est.CL_female$pr[1])), 
              H(est.CL_male$beta), (est.CL_male$sd.b)*(est.CL_male$pr[1]*(1-est.CL_male$pr[1])), NA, NA), 
            c(N_PL_female, est.PL_female$Se.Nhat, N_PL_male, est.PL_male$Se.Nhat, N_PL_tot, N_PL_tot.SE), 
            c(N_WPL_female, est.WPL_female$Se.Nhat, N_WPL_male, est.WPL_male$Se.Nhat, N_WPL_tot, N_WPL_tot.SE), 
            c(N_CL_female, est.CL_female$sd.Nhat, N_CL_male, est.CL_male$sd.Nhat, N_CL_tot, N_CL_tot.SE))

N_ests <- round(N_ests, digits = 3)
rownames(N_ests) <- c("PL", "WPL", "CL", "PL", "WPL", "CL")
colnames(N_ests) <- c("N_hat - female", "S.E.(N_hat)", "N_hat - male", "S.E.(N_hat)", "N_hat - combined", "S.E.(N_hat)")
N_ests

# Use full data and fit linear logistic regression (M_0/M_h) models.

# These results should be the same as Table 3.

# Partial likelihood approach:

# Construct PL weights to feed into glm().

R <- y - 1
h <- tau-t1
y.p <- R/h
y.p[is.na(y.p)] <- 0

est.PL_0 <- glm(y.p ~ 1, weights = h, family = binomial)
est.PL_const <- VarNhat.glm(est.PL_0, tau, y = y)

est.PL_1 <- glm(y.p ~ x.obs, weights = h, family = binomial)
est.PL_Mh <- VarNhat.glm(est.PL_1, tau, y = y)   

# Weighed partial likelihood approach:

# Construct WPL weights to feed into glm.

m.tilde.star <- tau-(tau + 1)/(y + 1)
h.wpl <- m.tilde.star
y.wpl <- R/h.wpl
y.wpl[is.na(y.wpl)] <- 0

est.WPL_0 <- glm(y.wpl ~ 1, weights = h.wpl, family = binomial)
est.WPL_const <- VarNhat.glm(est.WPL_0, tau, y = y)

est.WPL_1 <- glm(y.wpl ~ x.obs, weights = h.wpl, family = binomial)
est.WPL_Mh <- VarNhat.glm(est.WPL_1, tau, y = y) 

# Conditional likelihood approach.

Z <- matrix(c(rep(1, length(x.obs))), nrow = length(x.obs), ncol = 1)
est.CL_const <- par.mle(betain = 1, Z, tau, y = y)

Z <- matrix(c(rep(1,length(x.obs)), x.obs), nrow = length(x.obs), ncol = 2)
est.CL <- par.mle(betain = rep(1, ncol(Z)), Z, tau, y = y)

# Chao's lower-bound estimators and log-linear models (uses the Rcapture package).  

res_tot <- closedp(cap.hist)

# Combine results and display them.

N_ests <- matrix(round(as.numeric(rbind(rbind(est.PL_const[-1]), 
                                        rbind(est.WPL_const[-1]),
                                      c(est.CL_const$N.hat, est.CL_const$sd.Nhat),
                                      c(est.PL_Mh$Nhat, est.PL_Mh$Se.Nhat),
                                      c(est.WPL_Mh$Nhat, est.WPL_Mh$Se.Nhat),
                                      c(est.CL$N.hat,est.CL$sd.Nhat),
                                      c(res_tot$results[3, ][1], res_tot$results[3, ][2]),
                                      c(res_tot$results[4, ][1], res_tot$results[4, ][2]))),
                       digits = 3), ncol = 2)

rownames(N_ests) <- c("PL_0", "WPL_0", "CL_0", "PL_h", "WPL_h", "CL_h", "Chao", "Log-linear")

colnames(N_ests) <- c("N_hat", "S.E.(N_hat)")

round(N_ests, digits = 1)

#--------------------------------------------------------------------------------------------------------

rm(list = ls())

# Make sure all functions stored and loaded first in your R-console.

source("BiMJ_Functions_Hwang_et_al.R")

# Example 2: Variable selection using GLMNET and weighed partial likelihood.

# Results should be similar to Table 4.

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

dt2 <- dt[-which(dt$ofp == 0), ]

y <- dt2$ofp
n <- length(y)

X <- cbind(rep(1, n), dt2[, -1])
colnames(X)[1] <- "(Intercept)"

# Fit models and apply model selection (AIC and GLMNET).

# Fit a (full) log-linear positive-Poisson model using vglm().

dat1 <- data.frame(cbind(y, X))

mod1 <- vglm(y ~ emer + hosp + numchron + adldiff + age + black + gender + 
               married + school + faminc + employed + privins + medicaid, 
             pospoisson, data = dat1)

# Weighed partial likelihood approach:

# Construct WPL weights to feed into glm() and glmnet().

# Stepwise variable selection using AIC.

t.tilde.star <- y/(y + 1)
y.tilde <- y-1

dat2 <- data.frame(cbind(y.tilde, X))

offs <- log(t.tilde.star)

mod2 <- glm(y.tilde ~ emer + hosp + numchron + adldiff + age + black + 
              gender + married + school + faminc + employed + privins + medicaid,
            offset = offs, family = poisson, data = dat2)

AIC.glm2 <- stepAIC(mod2, trace = FALSE)$coefficients

# Elastic net variable selection using glmnet().

mod3 <- glmnet(as.matrix(X), y.tilde, family = "poisson", offset = log(t.tilde.star))

s <- cv.glmnet(as.matrix(X), y.tilde, family = "poisson", offset = log(t.tilde.star))$lambda.min

tmp_coeffs <- coef(mod3, s = s)

mod3.coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

# Display results for all three models.

coef(mod1)

AIC.glm2 

# Results below will be slightly different from the last column of Table 4 because glmnet() uses 
# cross-validation to select lambda, so the data is randomly split and will consist of 
# different training/test sets for each fit.

t(mod3.coef)

# K-fold cross validation.

K <- 5 

set.seed(2018) 

flds <- createFolds(y, k = K, list = TRUE, returnTrain = FALSE) 
names(flds)[1:K] <- rep("test", K) 

res.tr <- c() 
res.te <- c() 

for(kk in 1:K){
  set.seed(kk + 2018) 
  
  test <- flds[kk]$test   
  train <- (1:n)[-test] 
  
  # Construct training/test data.
  
  n.tr <- length(train) 
  n.te <- length(test) 
  
  y.tr <- y[train] 
  y.te <- y[test] 
  
  X.tr <- X[train, ] 
  X.te <- as.data.frame(X[test, ]) 
  
  y.tilde.tr <- y.tilde[train] 
  y.tilde.te <- y.tilde[test] 
  
  t.tilde.star.tr <- t.tilde.star[train] 
  t.tilde.star.te <- t.tilde.star[test] 
  
  est.names <- c("MLE (VGAM)", "WPL (GLM)", "WPL (GLMNET)") 
  
  # MLE log-linear VGAM.
  
  dat1 <- data.frame(cbind(y.tr, X.tr)) 
  
  a1 <- vglm(y.tr ~ emer + hosp + numchron + adldiff + age + black + gender + 
               married + school + faminc + employed + privins + medicaid, 
             pospoisson, data = dat1) 
  
  dat2 <- data.frame(X.te) 
  names(dat2) <- names(coef(a1)) 
  
  preds.tr0 <- c(exp(coef(a1)%*%t(as.matrix(dat1[, -1])))) 
  preds.tr <- preds.tr0/(1 - exp(-preds.tr0)) 
  var.y.tr <- preds.tr - preds.tr^2*exp(-preds.tr0) 
  
  preds.te0 <- c(exp(coef(a1)%*%t(as.matrix(dat2)))) 
  preds.te <- preds.te0/(1 - exp(-preds.te0)) 
  var.y.te <- preds.te - preds.te^2*exp(-preds.te0) 
  
  VGAM.DS.tr <- sum((y.tr - preds.tr)^2/var.y.tr + log(var.y.tr))/n.tr 
  VGAM.SPR.tr <- sum((y.tr - preds.tr)^2/var.y.tr)/n.tr 
  
  VGAM.DS.te <- sum((y.te - preds.te)^2/var.y.te + log(var.y.te))/n.te 
  VGAM.SPR.te <- sum((y.te - preds.te)^2/var.y.te)/n.te 
  
  # WPL.
  
  # Saturated log-linear (positive-Poisson MLE) GLM.
  
  offs <- log(t.tilde.star.tr) 
  est <- glm(y.tilde.tr ~ emer + hosp + numchron + adldiff + age + black + gender + 
               married + school + faminc + employed + privins + medicaid, 
             offset = offs, family = poisson, data = dat1) 
  
  AIC.glm <- stepAIC(est, trace = FALSE, direction = c("backward"))$coefficients 
  var.sel <- names(AIC.glm) 
  dat2a <- X.tr[, which(names(coef(est))%in%var.sel)] 
  preds.tr <- c(exp(AIC.glm%*%t(as.matrix(dat2a)))) 
  var.y.tr <- preds.tr-preds.tr^2*exp(-preds.tr0) 
  
  dat2b <- X.te[, which(names(coef(est))%in%var.sel)] 
  preds.te <- c(exp(AIC.glm%*%t(as.matrix(dat2b)))) 
  var.y.te <- preds.te - preds.te^2*exp(-preds.te0) 
  
  WPL.DS.tr <- sum((y.tr - preds.tr)^2/var.y.tr + log(var.y.tr))/n.tr 
  WPL.SPR.tr <- sum((y.tr - preds.tr)^2/var.y.tr)/n.tr 
  
  WPL.DS.te <- sum((y.te - preds.te)^2/var.y.te + log(var.y.te))/n.te 
  WPL.SPR.te <- sum((y.te - preds.te)^2/var.y.te)/n.te 
  
  # High-dimensional (elastic-net penalty) log-linear model using glmnet.  
  
  X.tr1 <- as.matrix(X.tr) 
  X.te1 <- as.matrix(X.te) 
  
  offs <- log(t.tilde.star.tr) 
  offset0.tr <- rep(0, n.tr) 
  offset0.te <- rep(0, n.te) 
  
  est <- glmnet(X.tr1[, -1], y.tilde.tr, family = "poisson", offset = log(t.tilde.star.tr)) 
  
  preds.tr0 <- predict(est, X.tr1[, -1], s = cv.glmnet(X.tr1[, -1], y.tilde.tr, family = "poisson", offset = log(t.tilde.star.tr))$lambda.min, type = "response", newoffset = offset0.tr) 
  preds.te0 <- predict(est, X.te1[, -1], s = cv.glmnet(X.tr1[, -1], y.tilde.tr, family = "poisson", offset = log(t.tilde.star.tr))$lambda.min, type = "response", newoffset = offset0.te) 
  
  preds.tr <- preds.tr0/(1 - exp(-preds.tr0))   
  preds.te <- preds.te0/(1 - exp(-preds.te0)) 
  
  var.y.tr <- preds.tr - preds.tr^2*exp(-preds.tr0) 
  var.y.te <- preds.te - preds.te^2*exp(-preds.te0) 
  
  WPL.GLMNET.DS.tr <- sum((y.tr - preds.tr)^2/var.y.tr + log(var.y.tr))/n.tr 
  WPL.GLMNET.SPR.tr <- sum((y.tr - preds.tr)^2/var.y.tr)/n.tr 
  
  WPL.GLMNET.DS.te <- sum((y.te - preds.te)^2/var.y.te + log(var.y.te))/n.te 
  WPL.GLMNET.SPR.te <- sum((y.te - preds.te)^2/var.y.te)/n.te 
  
  # Combine all results.
  
  DS1 <- c(VGAM.DS.tr, WPL.DS.tr, WPL.GLMNET.DS.tr) 
  SPR1 <- c(VGAM.SPR.tr, WPL.SPR.tr, WPL.GLMNET.SPR.tr) 
  
  DS2 <- c(VGAM.DS.te, WPL.DS.te, WPL.GLMNET.DS.te) 
  SPR2 <- c(VGAM.SPR.te, WPL.SPR.te, WPL.GLMNET.SPR.te) 
  
  res1 <- rbind(SPR1, DS1) 
  res2 <- rbind(SPR2, DS2) 
  
  colnames(res1) <- colnames(res2) <- c("positive-Poisson MLE", "WPL-AIC", "WPL-Elastic Net") 
  rownames(res1) <- rownames(res2) <- c("SPR", "DS") 
  
  res.tr <- c(res.tr, list(res1)) 
  res.te <- c(res.te, list(res2)) 
}

aa  =  bb <- matrix(0, nrow = nrow(res1), ncol = ncol(res1)) 

for(kk in 1:K){
  aa <- aa + res.tr[[kk]] 
  bb <- bb + res.te[[kk]] 
}

# Print results (bottom of Table 4).

est <- round(rbind(rbind(aa[1, ], bb[1, ])/K, rbind(aa[2, ], bb[2, ])/K), digits = 3)

row.names(est) <- c("SPR (train)", "SPR (test)", "DS (train)", "DS (test)") 

est
