set.seed(2018)

weight <- function(x) {
  return (dnorm(x, sd = sqrt(3)) / dt(x, 3))
}

for (S in c(100, 10000)) {
  t = rt(S, 3)
  w = unlist(lapply(t, weight))
  
  png(sprintf("1_weight_%d.png", S), width = 640, height = 480)
  hist(w, breaks = 100)
  dev.off()
  
  png(sprintf("1_log_weight_%d.png", S), width = 640, height = 480)
  hist(log(w), breaks = 100)
  dev.off()
  
  first.order = mean(t * w) / mean(w)
  second.order = mean(t^2 * w) / mean(w)
  
  expectation = first.order
  variance = second.order - first.order^2
  
  print(sprintf("1: size = %d, expectation = %f, 
                variance = %f", S, expectation, variance))
}

weight <- function(x) {
  return (dt(x, 3) / dnorm(x, sd = sqrt(3)))
}

for (S in c(100, 10000)) {
  n = rnorm(S, sd = sqrt(3))
  w = unlist(lapply(n, weight))
  
  png(sprintf("2_weight_%d.png", S), width = 640, height = 480)
  hist(w, breaks = 100)
  dev.off()
  
  png(sprintf("2_log_weight_%d.png", S), width = 640, height = 480)
  hist(log(w), breaks = 100)
  dev.off()
  
  first.order = mean(n * w) / mean(w)
  second.order = mean(n^2 * w) / mean(w)
  
  expectation = first.order
  variance = second.order - first.order^2
  
  print(sprintf("2: size = %d, expectation = %f, 
                variance = %f", S, expectation, variance))
}