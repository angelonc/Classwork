## Brain Size and Mental Capacity
#  Tasks: plot, correlate IQ with brain size, separately for men and women

rm(list = ls())
library(ggplot2)
library(ppcor)

# DIRECTORY STUFF (tab delimited w. header)
dir = '~/Documents/Classwork/Statistics_711/HW1/brain.dat'
dat = read.table(dir, header = 1, sep = '\t')

# PLOT
ggplot(dat, aes(x = MRI_Count, y = FSIQ, color = Female)) +
    geom_point() +
    xlab('Brain Size (voxels)') +
    ylab('Full Scale Intelligence Quotient (Wexler)') +
    geom_smooth(method = lm, se=0, fullrange=1)

# Factor by gender
dat$Female = factor(dat$Female, labels = c("Male","Female"))

# INGREDIENTS
# Means
dat$IQ.Mean = mean(dat$FSIQ)
dat$BS.Mean = mean(dat$MRI_Count)

# Deviation
dat$IQ.Diff = dat$FSIQ - dat$IQ.Mean
dat$BS.Diff = dat$MRI_Count - dat$BS.Mean
    
# Variance (squared deviance)
dat$IQ.Var = dat$IQ.Diff^2
dat$BS.Var = dat$BS.Diff^2

# Sum products (sum of deviations)
(SP = sum(dat$IQ.Diff * dat$BS.Diff))
# Sum of squares
(SS.IQ = sum(dat$IQ.Var))
(SS.BS = sum(dat$BS.Var))

# Descriptive stats
d.stats = apply(dat[,2:7], 2, function(x) c(
    Mean = mean(x, na.rm=1),
    Variance = var(x, na.rm=1),
    SD = sd(x, na.rm=1),
    Minimum = min(x, na.rm=1),
    Maximum = max(x, na.rm=1),
    Median = median(x,na.rm=1),
    Valid = sum(!is.na(x)),
    Missing = sum(is.na(x))
    )
)


dat$IQThresh = factor(dat$FSIQ > 110, labels = c('Lo', 'Hi'))
# Correlation
(SP / sqrt(SS.IQ*SS.BS))
cor.test(dat$FSIQ,dat$MRI_Count)
summary(lm(dat$FSIQ~dat$MRI_Count))

# (by gender)
by(dat, dat$Female, function(x) cor.test(x$FSIQ, x$MRI_Count))

# Plot lo vs hi
ggplot(dat, aes(x=MRI_Count, y = FSIQ, color = dat$IQThresh)) +
    geom_point() +
    geom_smooth(method = lm, se = 0, fullrange = 0)
by(dat, dat$IQThresh, function(x) cor.test(x$FSIQ, x$MRI_Count))

## IN CLASS
# Assign output to structure "Fit"
summary(Fit<-lm(dat$FSIQ~dat$MRI_Count))
# Expand structure contents
str(Fit)
# Plot it
hist(Fit$residuals)
# Find confidence intervals of estimates
confint(Fit)
# Plot estimates vs actual IQ
plot(dat$FSIQ~Fit$fitted.values)







