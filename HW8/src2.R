set.seed(2019)

#################### Slides 23-25
log_fc_alpha = function(theta, alpha, beta) {
  if (alpha<0) return(-Inf)
  n = length(theta)
  (alpha-1)*sum(log(theta))-n*lbeta(alpha,beta)-5/2*(alpha+beta)
}
log_fc_beta = function(theta, alpha, beta) {
  if (beta<0) return(-Inf)
  n = length(theta)
  (beta-1)*sum(log(1-theta))-n*lbeta(alpha,beta)-5/2*(alpha+beta)
}


mcmc = function(n_sims, dat, inits, tune0) {
  n_groups = nrow(dat)
  alpha = inits$alpha
  beta = inits$beta
  # Recording structure
  theta_keep = matrix(NA, nrow=n_sims, ncol=n_groups)
  alpha_keep = rep(alpha, n_sims)
  beta_keep = rep(beta , n_sims)
  tune = tune0
  for (i in 1:n_sims) {
    # Sample thetas
    theta = with(dat, rbeta(length(y), alpha+y, beta+n-y))
    # Sample alpha
    alpha_prop = rnorm(1, alpha, tune$alpha)
    logr = log_fc_alpha(theta, alpha_prop, beta)-log_fc_alpha(theta, alpha, beta)
    alpha = ifelse(log(runif(1))<logr, alpha_prop, alpha)
    # Sample beta
    beta_prop = rnorm(1, beta, tune$beta)
    logr = log_fc_beta(theta, alpha, beta_prop)-log_fc_beta(theta, alpha, beta)
    beta = ifelse(log(runif(1))<logr, beta_prop, beta)
    # Record parameter values
    theta_keep[i,] = theta
    alpha_keep[i] = alpha
    beta_keep[i] = beta
    
    # Adapt
    if (i > 20) {
      tune$alpha = 2.4^2 * var(alpha_keep[1:i]) / 2
      tune$beta = 2.4^2 * var(beta_keep[1:i]) / 2
    }
  }
  
  print(tune)
  
  return(data.frame(iteration=1:n_sims,
                    parameter=rep(c("alpha","beta",paste("theta[",1:n_groups,"]",sep="")),each=n_sims),
                    value=c(alpha_keep,beta_keep,theta_keep)))
}

d = read.table("./rats_data.txt",header=T)
dat=data.frame(y=d$y, n=d$n)
inits = list(alpha=1, beta=1)
# Run the MCMC
r = mcmc(2000, dat=dat, inits=inits, tune0=list(alpha=1,beta=1))

