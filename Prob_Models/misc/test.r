MUx = 10;
SIGx = sqrt(8);

x = (1:1000) * .3
y = (1:1000) * .2

pyx <- matrix(,nrow = 1000, ncol = 1000)
px = dnorm(x,MUx,SIGx)
for (i in length(y)) {
pyx[,i] <- dnorm(y[i], x, (.5 * .2 * x))
}
