# Set wd
setwd('~/Desktop/')

# Import data
dat <- read.csv("rclassdata.csv")
str(dat)

# Attach a data frame
attach(dat)
detach(dat)

# Operating on data frames
dat$iqpreham = (dat$iq * dat$ham_d_pre)
boys = subset(dat, sex == 1)
girls = dat[dat$sex == 2,]
t.test(boys$iq,girls$iq, var.equal = TRUE)

# For help:
?read.csv

# Installing a package
install.packages("robust_base")
