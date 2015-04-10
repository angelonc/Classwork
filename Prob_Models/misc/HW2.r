library(flux)
library(quantmod)

MUx <- 10
SIGx <- sqrt(8)

x <- (1:1000) * .03
y <- (1:1000) * .02

# Get prior distribution
p.x <- dnorm(x,MUx,SIGx)

# Preallocate
p.yx <- matrix(0,nrow=length(y),ncol=length(x))
p.y <- vector(mode='numeric',length=length(y))
p.xy <- matrix(0,nrow=length(y),ncol=length(x))
MAP <- vector(mode='numeric',length=length(y))

# For each y value, get posterior, likelihood, normalization
for (i in 1:length(y)) {
    p.yx[i,] <- dnorm(y[i], x, (0.5+0.2*x))
    p.y[i] <- auc(x,p.yx[i,]*p.x)
    p.xy[i,] <- p.yx[i,]*p.x/p.y[i]
    
    MAP[i] <- which(diff(sign(diff(p.xy[i,])))==-2)+1
}

plot(y, x[MAP])
points(x,x)
