# Code for samplers use in Beer sodium examples and the summary plots

beers <- read.table("beer_data.txt", header=T)
sodium.mat <- matrix(beers$Sodium,ncol=8,byrow=T)
sodium <- apply(sodium.mat,1,mean)

  
nsamp <- 5000
tau <- (0:4000)/100
sigma <- 0.846140/sqrt(8)
Vmu <- (sigma^2+tau^2)/6
muhat <- mean(sodium)

logpost <- tau*0

for (i in 1:length(tau)) {
  logpost[i] <- 0.5*log(Vmu[i]) -3*log(sigma^2+tau[i]^2) - sum((sodium-muhat)^2)/(sigma^2+tau[i]^2)/2
  }
post <- exp(logpost)

tausim <- sample(tau,nsamp,T,post)
Vmuk <- (sigma^2 + tausim^2)/6
musim <- rnorm(nsamp,muhat,sqrt(Vmuk))
thetasim <- matrix(0,ncol=6,nrow=nsamp)
for(i in 1:6) {
  vari <- sigma^2*tausim^2/(sigma^2 + tausim^2)
  meani <- (sodium[i]/sigma^2 + musim/tausim^2)*vari
  thetasim[,i] <- rnorm(nsamp,meani,sqrt(vari))
}

sodium.res <-sodium.mat
yrep <- array(0,c(nsamp,6,8))
r <- matrix(0, nrow=6, ncol=8)
corr <- rep(0,nsamp)
corrdat <- rep(0,nsamp)
vr <- rep(0,nsamp)
vrdat <- rep(0,nsamp)
for(i in 1:nsamp) {
  for(j in 1:6) {
    yrep[i,j,] <- rnorm(8,thetasim[i,j],sigma)
    r[j,] <- yrep[i,j,] - thetasim[i,j]
    sodium.res[j,] <- sodium.mat[j,] - thetasim[i,j]
  }
  corr[i] <- cor(sort(as.vector(r)),qnorm((1:48)/49))
  corrdat[i] <- cor(sort(as.vector(sodium.res)),qnorm((1:48)/49))
  vari <- apply(yrep[i,,],1,var)
  vr[i] <- max(vari)/min(vari)
  varidat <- apply(sodium.res,1,var)
  vrdat[i] <- max(varidat)/min(varidat)
}

sum(corr >= corrdat)/length(corr)
sum(vr >= vrdat)/length(vr)

fill.lm <- lm(Sodium ~ Beer, data=beers)
corrobs <- cor(sort(residuals(fill.lm)),qnorm((1:48)/49))

postscript("../beerdata.eps",horiz=F,width=6,height=4)
par(mfrow=c(1,1), pty="m", mar=c(4,4,2,1)+0.1)

plot(beers$Beer,beers$Sodium,pch=16, xlab="Beer", ylab="Sodium")
points(1:6,sodium,pch=16,col=2)
dev.off()

postscript("../beerhyper.eps",horiz=F,width=8,height=5)
par(mfrow=c(1,2), pty="m", mar=c(4,4,2,1)+0.1)

hist(tausim,prob=T,xlim=c(0,30), xlab=expression(tau),main="")
lines(tau,post/mean(post)/40,col=2)

hist(musim,prob=T, xlab=expression(mu),xlim=c(10,30),main="")

dev.off()

postscript("../beertheta.eps",horiz=F,width=8,height=5)
par(mfrow=c(2,3), pty="m", mar=c(4,4,2,1)+0.1)

hist(thetasim[,1],prob=T, xlab=expression(theta[1]),
  main="Beer 1")
points(sodium[1],0,pch=16,col=2)
points(mean(thetasim[,1]),0,pch=16,col=4)
hist(thetasim[,2],prob=T, xlab=expression(theta[2]),
  main="Beer 2")
points(sodium[2],0,pch=16,col=2)
points(mean(thetasim[,2]),0,pch=16,col=4)
hist(thetasim[,3],prob=T, xlab=expression(theta[3]),
  main="Beer 3")
points(sodium[3],0,pch=16,col=2)
points(mean(thetasim[,3]),0,pch=16,col=4)
hist(thetasim[,4],prob=T, xlab=expression(theta[4]),
  main="Beer 4")
points(sodium[4],0,pch=16,col=2)
points(mean(thetasim[,4]),0,pch=16,col=4)
hist(thetasim[,5],prob=T, xlab=expression(theta[5]),
  main="Beer 5")
points(sodium[5],0,pch=16,col=2)
points(mean(thetasim[,5]),0,pch=16,col=4)
hist(thetasim[,6],prob=T, xlab=expression(theta[6]),
  main="Beer 6")
points(sodium[6],0,pch=16,col=2)
points(mean(thetasim[,6]),0,pch=16,col=4)
dev.off()


postscript("../beercheck.eps",horiz=F,width=8,height=3.5)
par(mfrow=c(1,2), pty="m", mar=c(4,4,2,1)+0.1)

hist(corr,prob=T, xlab="Correlation",
  main="")
abline(v=corrobs,col=2)
hist(vr,prob=T, xlab="Variance Ratio",
  main="")
abline(v=vrdat,col=2)
dev.off()
