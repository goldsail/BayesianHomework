plot.increment = 0.0005
dat = read.csv("data.csv", header=T)

logit <- function(x) {
  return (log(x / (1 - x)))
}

log.odds = logit(dat$treated.deaths / dat$treated.total) - logit(dat$control.deaths / dat$control.total)
std.err = sqrt((dat$treated.deaths)^(-1) + 
                 (dat$treated.total - dat$treated.deaths)^(-1) + 
                 (dat$control.deaths)^(-1) +
                 (dat$control.total - dat$control.deaths)^(-1))

dat = cbind(dat, log.odds)
dat = cbind(dat, std.err)

J = nrow(dat)

# Question 5.15a

tau.posterior <- function(tau, dat) {
  J = nrow(dat)
  
  A = 0
  for (j in 1:J) {
    A = A + 1 / (dat$std.err[j]^2 + tau^2)
  }
  
  S1 = 0
  for (j in 1:J) {
    S1 = S1 + (dat$log.odds[j]) / (dat$std.err[j]^2 + tau^2)
  }
  
  S2 = 0
  for (j in 1:J) {
    S2 = S2 + (dat$log.odds[j]^2) / (dat$std.err[j]^2 + tau^2)
  }
  
  C = S2 - S1^2 / A
  
  P = 0
  for (j in 1:J) {
    P = P + (-1/2) * log(dat$std.err[j]^2 + tau^2)
  }
  P = exp(P)
  
  return (P * A^(-1/2) * exp(-C/2))
  
}

taus = seq(0 + plot.increment, 1, by = plot.increment)
taus.posterior = unlist(lapply(taus, tau.posterior, dat))
taus.posterior = taus.posterior / sum(taus.posterior) / plot.increment # rescale to 1

png("tau.posterior.png", width = 800, height = 600)
plot(taus, taus.posterior)
dev.off()

# Question 5.15b

theta.posterior.mean <- function(tau, dat) {
  J = nrow(dat)
  
  S0 = 0
  for (j in 1:J) {
    S0 = S0 + 1 / (dat$std.err[j]^2 + tau^2)
  }
  
  S1 = 0
  for (j in 1:J) {
    S1 = S1 + (dat$log.odds[j]) / (dat$std.err[j]^2 + tau^2)
  }
  
  ret = 1:12
  for (j in 1:J) {
    A = dat$log.odds[j] / dat$std.err[j]^2 + 1 / tau^2 * S1 / S0
    B = 1 / dat$std.err[j]^2 + 1 / tau^2
    ret[j] = A / B
  }
  
  return (ret)
}

theta.posterior.var <- function(tau, dat) {
  J = nrow(dat)
  
  S0 = 0
  for (j in 1:J) {
    S0 = S0 + 1 / (dat$std.err[j]^2 + tau^2)
  }
  
  ret = 1:12
  for (j in 1:J) {
    A = 1 / tau^2
    B = 1 / dat$std.err[j]^2 + 1 / tau^2
    ret[j] = 1 / B + (A / B)^2 / S0
  }
  
  return (ret)
}

theta.posterior.means = lapply(taus, theta.posterior.mean, dat)
theta.posterior.means = matrix(unlist(theta.posterior.means), ncol=J, byrow=T)

theta.posterior.vars = lapply(taus, theta.posterior.var, dat)
theta.posterior.vars = matrix(unlist(theta.posterior.vars), ncol=J, byrow=T)

png("theta.posterior.means.png", width = 800, height = 600)
matplot(taus, theta.posterior.means, type="l", lty=1, main="Posterior means of theta")
dev.off()

png("theta.posterior.vars.png", width = 800, height = 600)
matplot(taus, theta.posterior.vars, type="l", lty=1, main="Posterior variances of theta")
dev.off()

# Question 5.15c

index.MLE = which.max(taus.posterior) # use MLE of tau for a crude estimate of theta_j
tau.MLE = taus[index.MLE]

# Question 5.15d

p.mu.tau.c.dat <- function(mu, tau, dat) {
  J = nrow(dat)
  
  A = 0
  for (j in 1:J) {
    A = A + 1 / (dat$std.err[j]^2 + tau^2)
  }
  
  S1 = 0
  for (j in 1:J) {
    S1 = S1 + (dat$log.odds[j]) / (dat$std.err[j]^2 + tau^2)
  }
  
  S2 = 0
  for (j in 1:J) {
    S2 = S2 + (dat$log.odds[j]^2) / (dat$std.err[j]^2 + tau^2)
  }
  
  C = S2 - S1^2 / A
  
  S3 = 0
  for (j in 1:J) {
    S3 = S3 + (dat$log.odds[j]) / (dat$std.err[j]^2 + tau^2)
  }
  
  B = S3 / A
  
  P = 0
  for (j in 1:J) {
    P = P + (-1/2) * log(dat$std.err[j]^2 + tau^2)
  }
  P = exp(P)
  
  return (P * exp(-A/2 * (mu - B)^2) * exp(-C/2))
}

p.theta.c.mu.tau.dat <- function(j, theta.j, mu, tau, dat) {
  J = nrow(dat)
  
  A = 1 / dat$std.err[j]^2
  B = 1 / tau^2
  
  theta.j.hat = (A * dat$log.odds[j] + B * mu) / (A + B)
  V.j = 1 / (A + B)
  
  return (V.j^(-1/2) * exp(-1/2 * (theta.j - theta.j.hat)^2 / V.j))
}

# Question 5.15d

T = 1000
samples = matrix(NA, T, J)

for (k in 1:1000) {
  U = runif(1)
  temp = 1:length(taus)
  temp[1] = taus.posterior[1]
  V = 0
  for (j in 2:length(taus)) {
    temp[j] = temp[j - 1] + taus.posterior[j]
    if (temp[j] > U) {
      tau.sampled = taus[j] # sample tau
      break;
    }
  }
  
  J = nrow(dat)
  A = 0
  for (j in 1:J) {
    A = A + 1 / (dat$std.err[j]^2 + tau.sampled^2)
  }
  S3 = 0
  for (j in 1:J) {
    S3 = S3 + (dat$log.odds[j]) / (dat$std.err[j]^2 + tau.sampled^2)
  }
  B = S3 / A
  mu.sampled = rnorm(1, mean = B, sd = 1 / sqrt(A)) # sample mu
  
  J = nrow(dat)
  for (j in 1:J) {
    C = 1 / dat$std.err[j]^2
    D = 1 / tau.sampled^2
    
    theta.j.hat = (C * dat$log.odds[j] + D * mu.sampled) / (C + D)
    V.j = 1 / (C + D)
    theta.j.sampled = rnorm(1, theta.j.hat, sqrt(V.j)) # sample theta_j
    
    samples[k, j] = theta.j.sampled
  }
}

png("theta.simulations.png", width = 800, height = 600)
par(mfrow=c(4,6))
for (j in 1:J) {
  hist(samples[, j], main="Simulations of theta", pch = j)
}
dev.off()