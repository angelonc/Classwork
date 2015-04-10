#### PSYC 739 HW#5
#### Prof Alan Stocker, Spring 2015
#### Code by: Welton Chang
install.packages("pracma")
library(pracma)

### p(x|m) = p(x|m)p(m) / auc (p(x|m)p(m)) dθ 
# p(x|m)*p(m)

# A.1.2 #

#### set value for x
x<- linspace(1, 100, 1000)

#### draw sensory signal m of test line segments, mu = 70, sd = 2.5 
m <- rnorm(1000, 70, 2.5) 
hist(m)

#### Draw prior distribution

prior <- dnorm(x, 50, 5)

### Pre allocate

m1 <- vector(mode='numeric', length(x))
posterior <- vector(mode='numeric', length(x))
posteriorm <-matrix(0, nrow=1000, ncol=1000)
BLS <- vector(mode='numeric', length(x))

#### for loop
for (i in 1:length(m)) {
  m1 <- dnorm(m[i], x, 2.5)
  posterior <- (prior*m1) / auc(x, prior*m1)  
  BLS[i] <- auc(x, posterior*x)
}

hist(BLS, xlab="BLS values", main="BLS")
hist(BLS, freq=FALSE, xlab="BLS values", main="BLS")


#### Find BLS

s <- dnorm(x, 70, 2.5)
hist(BLS, freq=FALSE)
lines(x, s, col="red")
lines(x, prior, col="blue", lwd=3)
lines(x, -1*((prior*s)/auc(x, prior*s)), col="green", lwd=3)


#### A.1.4

### Pre allocate

xr <- linspace(55, 85, 13)
response <- matrix(0, nrow=length(xr), ncol=100)

for(i in 1:length(xr)) {
  mtest <- rnorm(100, 70, 2.5)
  mref <- rnorm(100, xr[i], 1.5)
  response[i,] <- mref>mtest   
}

probreftest<- apply(response,1,mean)

plot(xr,probreftest, col="blue", main="Psychometric curve 2AFC task A.2", xlab="values", ylab="P(r>t)")
lines(xr,probreftest, col="red", lwd=3)

#######
# A.2.2

posteriorPDF <-matrix(0, nrow=1000, ncol=1000)
probrefgtest <- vector(mode='numeric', length(xr))

for(i in 1:length(xr)) {
  mrefPDF <- dnorm(x, xr[i], 1.5)
  mtestCDF <- pnorm(xr[i], 70, 2.5)
  mrefposteriorPDF <- -1*((mrefPDF*prior)/auc(x, prior*mrefPDF))
  probrefgtest[i] <- -1*auc(x, mrefposteriorPDF)*mtestCDF
}

for(i in 1:length(x)) {
  mrefPDF <- dnorm(x, x[i], 1.5)
  mtestCDF <- pnorm(x[i], 70, 2.5)
  mrefposteriorPDF <- -1*((mrefPDF*prior)/auc(x, prior*mrefPDF))
  probrefgtest[i] <- -1*auc(x, mrefposteriorPDF)*mtestCDF
}

#### Plot results

plot(x, probrefgtest, col="blue", main="Precise values of p(xhat|x)", xlab="Values", ylab="Probability")

