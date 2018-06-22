library("arm")
library("BMS")
library("bayess")
library("BayesFactor")

# Q1

dat = datasets::longley

y = dat$GNP.deflator
X = dat$Population

mod = arm::bayesglm(y ~ X, 
                    prior.mean=0,
                    prior.scale=Inf,
                    prior.df=Inf)

summary(mod)

n = length(y)
k = 2
s2 = sum(mod$residuals ^ 2)/(n-k)

posterior.mean.sigmasq = (n-k) * s2 / (n-k-2)
posterior.mean.sigmasq

# Q2 

mod = arm::bayesglm(y ~ X, 
                    prior.mean=0,
                    prior.scale=1,
                    prior.df=Inf,
                    prior.mean.for.intercept=0,
                    prior.scale.for.intercept=1,
                    prior.df.for.intercept=Inf)

summary(mod)
posterior.mean.beta = mod$coefficients
posterior.mean.beta

# Q3

mod.full = BayesReg(dat$GNP.deflator, dat[, -1], g=nrow(dat))
mod = BayesReg(dat$GNP.deflator, dat[, c(-1, -3)], g=nrow(dat))

mod.full
mod

bf.full = lmBF(GNP.deflator ~ 
                 GNP + Unemployed + Armed.Forces + Population + Year + Employed, 
               data = dat)
bf = lmBF(GNP.deflator ~ 
            GNP + Armed.Forces + Population + Year + Employed, 
          data = dat)

bf.full / bf
