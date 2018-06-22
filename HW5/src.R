# Question 6.1
data.observed = c(28, 18, 12, 8, 7, 1, -1, -3)
data.simulated = c(26, 19, 14, 10, 6, 2, -3, -9)

png("data.vs.png", width=800, height=600)
plot(data.observed, data.simulated, main="Simulated v.s. simulated order statistics")
dev.off()

mod = lm((data.simulated - data.observed) ~ data.observed)
summary(mod)


# Question 6.6
T0 = 3 # number of switches
N = 10000
T = 1:N
theta = rbeta(N, shape1=7+1, shape2=13+1)
for (j in 1:N) {
  count.switches = 0
  y.last = rbinom(1, size=1, prob=theta[j])
  count.zeros = 1 - y.last
  while (1) {
    y = rbinom(1, size=1, prob=theta[j])
    if (y != y.last) {
      count.switches = count.switches + 1
    }
    count.zeros = count.zeros + 1 - y.last
    if (count.zeros > 12.5) {
      break;
    }
    y.last = y
  }
  T[j] = count.switches
}
png("switches.histogram.png", width=800, height=600)
hist(T, breaks=seq(-0.5, max(T)+0.5))
dev.off()