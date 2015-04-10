##############################
### IN CLASS DEMO (9/9/14) ###
##############################

# Clear all
rm(list = ls())

x <- c(1,2,2,3,4) # <- is an assignment operator, c() concatenates
mean(x)
sd(x)

sum((x - mean(x))^2)
var(x)
median(x)
