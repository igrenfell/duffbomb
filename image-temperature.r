setwd("I:\\workspace\\restime")


resmat <- read.csv("ResTimeData.csv")

resmat <- read.csv("ParticleSpacing_Data_Correctloadings.csv")
resmat <- read.csv("ParticleSpacing_Data_samedensities.csv")

resmat <- resmat[1:50,]
restime2 <- resmat$Residence.Time..s.
Betavals1 <- resmat$beta1..Packing.Ratio...mm.3.mm.3.
Betavals2 <- resmat$beta2..Packing.Ratio...mm.3.mm.3.
Betavals3 <- resmat$beta3..Packing.Ratio...mm.3.mm.3.
lambdavals <- resmat$Lamda..Porosity...mm.3.mm.2.

loadingvec1 <- resmat$Fuel.Loading..kg.m.2.
loadingvec2 <- resmat$Fuel.Loading..g.mm.2.
loadingvec3 <- resmat$Fuel.Loading..tons.acre.

diam <- resmat$Diameter..mm.

massvec <- resmat$mass..g.
lrt <- log(restime2)
logbeta2 <- log(Betavals2)

squaremod.linear <- lm(restime2~ diam * poly(Betavals2, 2, raw=TRUE))
summary(squaremod.linear)

plot(Betavals2, restime2, type = "n")
betapredseq <- seq(min(Betavals2), max(Betavals2), length= 100)
diam.unique <- unique(diam)
logdiam.unique <- log(diam)
ndiams <- length(diam.unique)

prd <- Betavals2**2 / diam

for(i in 1:ndiams)
{
  predmat <- data.frame(Betavals2 = betapredseq, diam = rep(diam.unique[i], 100))
  predvec <- predict(squaremod.linear, predmat)
  lines((betapredseq), (predvec), col  = i)
  points((Betavals2[diam == diam.unique[i]]), (restime2[diam == diam.unique[i]]), col = i)
  
}

plot(squaremod.linear$residuals)

plot(prd, restime2)
abline(coef(simple.model)[1], coef(simple.model)[2])

transform.model <- lm(restime2~  diam+poly(prd, 3, raw=TRUE))
summary(transform.model)


plot(prd, restime2, type ="n")
prdpredseq <- seq(min(prd), max(prd), length= 100)

prd.sub <- prd[diam ==diam.unique[3]]
lrt.sub <- lrt[diam ==diam.unique[3]]
# transform.model <- lm(restime2~  diam+prd)
# summary(transform.model)

prd <- (diam)*((Betavals2^0.3333))
logdiam <- log(diam)
prd <- (Betavals2^2) / logdiam

dm <- data.frame(Betavals = Betavals2, logdiam = logdiam)
nlsmod <- nls(lrt ~ A*Betavals^k / logdiam + B * logdiam, start = list(A =1, B = 1, k = 2), data = dm)
summary(nlsmod)

simple.model <- lm(lrt ~ prd + loadingvec1 + logdiam)


#prd <- (diam)*((Betavals2^0.3333))
prdpredseq <- seq(min(prd), max(prd), length= 100)
loadpredseq <- seq(min(loadingvec1), max(loadingvec1), length= 100)

betapredseq <- seq(min(Betavals2), max(Betavals2), length= 100)

logmass <- log(massvec)
simple.model <- lm(lrt ~ loadingvec1 + logdiam)
summary(simple.model)
plot(loadingvec1, lrt, type  = "n", ylim = c(0, 6), xlab = "Loading (kg / m^2)")
#plot(Betavals2, lrt, type  = "n", ylim = c(0, 6))

for(i in 1:ndiams)
{
 # predmat <- data.frame(prd = prdpredseq, logdiam = log(rep((diam.unique[i]), 100)), loadingvec = rep(median(loadingvec), 100))
  predmat <- data.frame(loadingvec1 = loadpredseq, logdiam = log(rep((diam.unique[i]), 100)))
  
  #predmat <- data.frame(Betavals = betapredseq, logdiam = rep(log(diam.unique[i]), 100))
  
  predvec <- predict(simple.model, predmat)
  #predvec <- predict(nlsmod, predmat)
  lines((loadpredseq), (predvec), col  = i)
  points(loadingvec1[diam == (diam.unique[i])], lrt[diam == diam.unique[i]], col = i)
  
}

summary(simple.model$residuals)


models <- list()
plot(prd, lrt, type  = "n")

for(i in 1:ndiams)
  {
  # Subset data for the current diameter
  prd.sub <- prd[diam == diam.unique[i]]
  lrt.sub <- lrt[diam == diam.unique[i]]
  # Fit a polynomial regression (e.g., degree 2)
  #poly_model <- lm(lrt.sub ~ poly(prd.sub, 2))
  
  # Store the model in the list
  #models[[as.character(diam.unique[i])]] <- poly_model
  
  # Print a summary for inspection
  points(prd.sub, lrt.sub, col = i)
  prdpredseq <- seq(min(prd.sub), max(prd.sub), length= 100)
  
  predmat <- data.frame(prd = prdpredseq, diam = rep(diam.unique[i], 100))
  
  lines(prdpredseq, predict(fm1, predmat), col  =i)
}

dfac <- as.factor(diam)
fm1 <- lmer(lrt ~ poly(prd, 3) + ( poly(prd, 3) | diam))
summary(fm1)

squaremod <- lm(lrt ~ diam * poly(logbeta2, 2, raw=TRUE))

summary(squaremod)

predmat <- data.frame(prd = prdpredseq, dfac = rep(diam.unique[1], 100))
predvec <- predict(fm1, predmat)

plot(prdpredseq, predvec)

cubemod <- lm(lrt ~ diam * poly(logbeta2, 3, raw=TRUE))

summary(cubemod)

betapredseq <- seq(min(logbeta2), max(logbeta2), length= 100)
diam.unique <- unique(diam)
ndiams <- length(diam.unique)
plot(logbeta2, lrt, type = "n")

for(i in 1:ndiams)
{
  predmat <- data.frame(logbeta2 = betapredseq, diam = rep(diam.unique[i], 100))
  predvec <- predict(squaremod, predmat)
  lines((betapredseq), (predvec), col  = i)
  points((logbeta2[diam == diam.unique[i]]), (lrt[diam == diam.unique[i]]), col = i)
  
}

plot(logbeta2, lrt)

for(i in 1:ndiams)
{
  predmat <- data.frame(logbeta2 = betapredseq, diam = rep(diam.unique[i], 100))
  predvec <- predict(cubemod, predmat)
  lines(betapredseq, predvec, col  = i)
  points(logbeta2[diam == diam.unique[i]], lrt[diam == diam.unique[i]], col = i)
  
}

Diam <- resmat$Diam
Loading <- resmat$loading
Moisture <- resmat$Moisture

nlsmod <- nls(restime2 ~ A*Betavals + B*Diam^DiamExp + C*Loading + D*Moisture + K, start = list(K = -2.0, A = 571, B = 57000, C = 3.0, D = 0.1, DiamExp = 1.5))

summary(nlsmod)


rt2 <- read.csv("CopyofRysData.csv")

rtcomb <- read.csv("CombBurn_ResTimeData.csv")

restime2 <- rt2$Residence.Time..s.
Betavals <- rt2$beta..Packing.Ratio...mm.3.mm.3.
Diam <- rt2$Diameter..mm.
Loading <- rt2$Fuel.Loading..tons.acre.



nlsmod <- nls(restime2 ~ A*Betavals + B*Diam^DiamExp + C*Loading + K, start = list(K = -2.0, A = 571, B = 57000, C = 3.0, DiamExp = 1.5))

summary(nlsmod)

plot(log(Diam), log(restime))

logdiam <- log(Diam)
lrt <- log(restime2)
plot(Loading, lrt)

plot(Betavals, lrt)

logdiam <- log(Diam)
logloading <- log(Loading)

lmod <- lm(lrt ~ logloading  +logdiam + Betavals *logdiam )
summary(lmod)

plot(lmod$residuals)

betaseq <- seq(0, 1, by = .01)

a1 <- -1
a2 <- 0.5
betacrit <- 0.5
b1 <- -10
b2 <- 10
d <- 3
restimemod <- a1*exp(b1*(betaseq - betacrit))+a2*exp(b2*(betaseq - betacrit))+d
plot(betaseq, restimemod, type = "l")


a1 <- 0
b1 <- 1
c1 <- 2
e1 <- 0.1
f1 <- 1

restimemod2 <- a1+b1*(betaseq - betacrit)+c1*(betaseq - betacrit)^2 +e1*exp(f1*betaseq - betacrit)
plot(betaseq, restimemod2, type = "l")


plot(betaseq, restimemod, type = "l")


lfit <- predict(lmod)
plot(lrt, lfit)
abline(0, 1)

head(rtcomb)

diam_comb <- rtcomb$Diam..mm.
loading_comb <- rtcomb$loading..kg.m2.
rt_comb <- rtcomb$restime..s.
rt_beta <- rtcomb$Beta

logdiam_comb <- log(diam_comb)
logloading_comb <- log(loading_comb)
logrt_comb <- log(rt_comb)

lmod_comb <- lm(logrt_comb ~ logloading_comb * logdiam_comb + rt_beta*logdiam_comb)
summary(lmod_comb)


logrt_all <- c(lrt, logrt_comb)
logloading_all <- c(logloading, logloading_comb)
logdiam_all <- c(logdiam, logdiam_comb)
beta_all <- c(Betavals, rt_beta)
iscomb <- c(rep(0, length(lrt)), rep(1, length(logrt_comb)))

allmod <- lm(logrt_all ~ logloading_all*logdiam_all+  beta_all*logdiam_all)
summary(allmod)

allpred <- predict(allmod)
plot(logrt_all, allpred)

points(logrt_all[iscomb == 0], allpred[iscomb == 0], col = "red")

abline(0, 1)






