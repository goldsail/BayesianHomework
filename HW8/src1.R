set.seed(2018)
library(sos)

dat = cbind(
  c(83,92,92,46,67),
  c(117,109,114,104,87),
  c(101,93,92,86,67),
  c(105,119,116,102,116),
  c(79,97,103,79,92),
  c(57,92,104,77,100)
)

n = 5
N = 30

MH <- function(f, v, x0) {
  # f: the target pdf
  # v: variance of the proposal normal distribution
  # x0: initial point
  
  N = 5000
  x = 1:N
  x[1] = x0
  for (i in 2:N) {
    t = rnorm(1, x[i-1], v)
    r = f(t) / f(x[i-1])
    u = runif(1)
    if (r > u) {
      x[i] = t
    } else {
      x[i] = x[i-1]
    }
  }
  return(x[N])
}

# Separate Model

gibbs_separate <- function(theta_0, n_points, y) {
  theta = matrix(theta_0, nrow=n_points, ncol=2, byrow=T)
  n = length(y)
  for (i in 2:n_points) {
    theta[i, 1] = rnorm(1, mean(y), theta[i-1, 2] / n)
    theta[i, 2] = MH(function(x) {
        return(x^{-1-n/2} * exp(-1/(2*x) * ((n-1)*sd(y)^2 + n*(mean(y)-x)^2)))
      }, 5, 100)
  }
  return(theta)
}

n_points = 500
png("posterior_mean_separate.png", width = 640, height = 480)
theta_separate = gibbs_separate(100, n_points, y=dat[,6])
hist(theta_separate[,1], breaks=20, main="Histogram of posterior mean (separate)")
dev.off()

# Pooled Model

n_points = 500
png("posterior_mean_pooled.png", width = 640, height = 480)
theta_pooled = gibbs_separate(100, n_points, y=as.vector(dat))
hist(theta_pooled[,1], breaks=20, main="Histogram of posterior mean (pooled)")
dev.off()

h_pooled = matrix(0, nrow=n_points, ncol=10, byrow=T)
for (i in 1:n_points) {
  h_pooled[i,] = rnorm(10, mean=theta_pooled[i,1], sd=sqrt(theta_pooled[i,2]))
}
png("posterior_predictive_pooled.png", width = 640, height = 480)
hist(as.vector(h_pooled), breaks=20, main="Histogram of posterior predictive (pooled)")
dev.off()
png("posterior_predictive_mean_pooled.png", width = 640, height = 480)
hist(rowMeans(h_pooled), breaks=20, main="Histogram of posterior predictive (pooled)")
dev.off()

# Hierarchical Model

gibbs_hierarchical <- function(theta_0, n_points, y) {
  theta = matrix(theta_0, nrow=n_points, ncol=9, byrow=T)
  n = length(y)
  n_j = 5
  for (i in 2:n_points) {
    
    for (j in 1:6) {
      theta[i,j] = rnorm(1, 
        mean=(theta[i-1,7]*1/theta[i-1,9] + mean(y[,j])*n_j/theta[i-1,8]) / (1/theta[i-1,9] + n_j/theta[i-1,8]),
        sd=sqrt(1/(1/theta[i-1,9] + n_j/theta[i-1,8])))
    }
    
    theta[i,7] = rnorm(1, mean=mean(theta[i,1:6]), sd=sqrt(theta[i-1, 9]/6))
    
    temp = 0
    for (k in 1:6) {
      for (u in 1:n_j) {
        temp = temp + (y[u,k] - theta[i, k])^2
      }
    }
    theta[i,8] = 1 / rgamma(1, 6*n_j/2, (temp/(6*n_j)/2))
    
    temp = 0
    for (k in 1:6) {
      temp = temp + (theta[i, k] - theta[i, 7])^2
    }
    theta[i,9] = 1 / rgamma(1, (6-1)/2, (temp/(6-1)/2))
  }
  return(theta)
}

n_points = 50000
theta_hierarchical = gibbs_hierarchical(100, n_points, y=dat)
png("posterior_mean_hierarchical.png", width = 640, height = 480)
hist(theta_hierarchical[,6], breaks=50, main="Histogram of posterior mean (hierarchical, for 6th machine)")
dev.off()

h_hierarchical = matrix(0, nrow=n_points, ncol=10, byrow=T)
for (i in 1:n_points) {
  h_hierarchical[i,] = rnorm(10, 
                             mean=rnorm(1, theta_hierarchical[i,7], sqrt(theta_hierarchical[i,9])), 
                             sd=sqrt(theta_hierarchical[i,8]))
}
png("posterior_predictive_hierarchical.png", width = 640, height = 480)
hist(as.vector(h_hierarchical), breaks=50, main="Histogram of posterior predictive (hierarchical)")
dev.off()
png("posterior_predictive_mean_hierarchical.png", width = 640, height = 480)
hist(rowMeans(h_hierarchical), breaks=50, main="Histogram of posterior predictive (hierarchical)")
dev.off()