######################
## IN CLASS R STUFF ##
######################

# Data
X <- c(8,8,8,8,8,8,8,19,8,8,8)
Y <- c(6.58,5.76,7.71,8.84,8.47,7.04,5.25,12.5,5.56,7.91,6.89)

# Stats
sums = c(sum(X),sum(Y))
means = c(mean(X), mean(Y))
stds = c(sd(X),sd(Y))
r = cor.test(X,Y)
r2 =cor(X,Y)^2
LM = summary(lm(Y~X))

# Plot
plot(X,Y)

## Blah
install.packages('ppcor', depend=1)
library('ppcor')

x1 = c(4,4,7,7,10,10)
x2 = c(1,2,2,4,3,6)
y = c(14,23,30,50,39,67)

cor(data.frame(x1,y,x2))

cor.test(x1,x2)                # r = .75
cor.test(x1,y)                 # r = .81
cor.test(x2,y)                 # r = .99

pcor.test(x1,y,x2)             # r = .62
pcor.test(x2,y,x1)             # r = .98

spcor.test(x1,y,x2)            # r = .41
spcor.test(x2,y,x1)            # r = .65

# Effect of x2 on x1
fit <- lm(x1 ~ x2)
cor(fit$residuals, y)          # r = .09
cor(fit$residuals, y)^2        # r2 = .009

# Effect of x1 on x2 
fit1 <- lm(x2 ~ x1)
cor(fit1$residuals, y)         # r = .58
cor(fit1$residuals, y)^2       # r2 = .34

fit2 <- lm(x1 ~ y)
cor(fit2$residuals, x2)
cor(fit2$residuals, y) # Literally removed effect of y from x1

spcor(data.frame(x1,y,x2))
pcor(data.frame(x1,y,x2))

##########

rm(list = ls())

setwd('~/Documents/Classwork/Statistics_711/')

install.packages('car',depend=TRUE)

library(car)
library(ppcor)

x = read.table('Classdemo.dat', header=1)
head(x) # looks at first pary of data
x[1:6, ] # same thing

Fit1 = lm(y ~ x1 + x2 + x3 + x4, data = x)
summary(Fit1) # testing whether all of these variables differ from 0

###########################
### MULTIPLE REGRESSION ###
###########################

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.008331   0.018589   0.448    0.654     # 65% chance that you would see an abs value of t = .448
#x1          0.955822   0.111985   8.535  < 2e-16 ***
#x2          0.877780   0.109472   8.018 1.54e-15 ***
#x3          0.936303   0.111937   8.365  < 2e-16 ***
#x4          0.996143   0.080694  12.345  < 2e-16 ***

# Probability of observing a t-statistic of at least t-value (65%)
# .008 * +/-(.02*1.96) = CI

# lm is always semi-partial

# X1 Summary:
# .96 change in y for each unit x1 holding all things constant
# .11 stddev of sample dist
# If null is true, it is very unlikely that we would find a t-value as large as we did... x1 is related to outcome (holding everything else constant)

#Residual standard error: 1.001 on 2896 degrees of freedom
#Multiple R-squared:   0.05,	Adjusted R-squared:  0.04869 
#F-statistic: 38.11 on 4 and 2896 DF,  p-value: < 2.2e-16

# Residual Standard Error: the standard deviation of the errors of all the variables around the regress line
# Multiple r-squared: correlation of fitted values vs observed values
# Adjusted r-squared: penalizes you for the number of parameters in the model
# df = 2901 w. 5 estimates = 2896 df
# F-statistic: 38.11 on 4 (variables) and 2896 df

spcor.test(x$y, x$x1, x[ , c('x2','x3','x4')])

Fit2 = lm(x1 ~ x2 + x3 + x4, data=x)
cor( Fit2$residuals, x$y)

cor(x$y, Fit1$fitted.values)^2 # correlation between actual and predicted values (same as multiple r-squared)

# All regression is doing is trying to come up with weights so that fitted values are maximally correlated with the actual values

hist(Fit1$residuals, las=1, xlab="Residual")
plot(Fit1$residuals ~ Fit1$fitted.values, pch=19) # always plot residuals

######################################
### VIF: variance inflation factor ###
######################################

vif(Fit1)

# With shared variance c you don't know who to give credit to
# VIF are the degree to which standard errors have increased due to multicolinearity with other predictors
# If predictors were completely independent, VIFs would be 1
# How big is big in VIFs?
# Worrysome around 8 9 or 10

Fit3 = lm(x1 ~ x2 + x3 + x4, data=x)
summary(Fit3) 
tol = 1 - summary(Fit3)$r.squared # tolerance (how resistant x1 is to variance of other variables)
# x2 3 and 4 in combination predict about 47% of variability in x1 - x1 has about 53% variance independent from the other variables
1/tol
1/vif(Fit1)

# Because of the multicolinearity of the other variables, the standard error of x1 will increase by a factor of 1.899 -- inverse to find r-squared


#############
### PLOTS ###
#############
residualPlots(Fit1)[-5]
# We want splines to be flat -- no relationship betweenn predictors and residuals
avPlots(Fit1, id.n=2)
influencePlot(Fit1, id.n=3)
# Studentized residuals: put residuals on a t-scale (look for outliers)
# Hat-Values: values on the far right could potentially change the slope and intercept - go back and look at the labeled observations -- double-check source data to determine why it is a weird observation (often data-entry error) - can drop it because of its influence on the fit (maybe a new model is more appropriate though)
# Leverage: outliers on the y-values... the farther an observation is from the center the more influence it has on the data (see-saw metaphor)
#  - outlier can change slope and intercept of fitted value (1 value can change the slope and intercept)
# Influential observations: outliers on the x-values
#  - can be an outlier but still not be influential
# Finds influence by dropping one data point at a time and rerunning the analysis to see how the estimates are changed for each variable removed
ncvTest(Fit1) # statistical test for homoskedasticity

confint(Fit1)
# Better than p-value
# - tells us stability of estimates
# - tells us if our sample estimate is significant


# HOMEWORK
# Analyze regression data
# WRite up results
# Think about whether data satisfy assumptions of linear model
# Be prepared to talk about our healthcare related behaviors
# Consequences of what you observe in your data when diagnosing your model



############## Tues Nov 04th ##################
# One factor, fixed effects ANOVA
rm (list = ls())
library(ggplot2)

Temperature = c(1,1,1,1,2,2,2,2,3,3,3,3)
Correct =     c(2,4,5,1,9,8,5,10,6,3,4,3)

X = data.frame(Temperature, Correct)

X$Temperature = factor(X$Temperature, levels = 1:3, labels = c('50 deg.', '70 deg.', '90 deg.'))

aggregate(Correct ~ Temperature, X, FUN = function(x) c(Mean = mean(x), SD = sd(x), Median = median(x)))

ggplot(X, aes(x=Temperature, y = Correct, fill = Temperature)) +
       geom_boxplot() +
       xlab('Temperature') +
       ylab('Correct Responses')

qf(.95, 2, 9)

# ANOVA
summary(AOV1 <- aov(Correct ~ Temperature, data=X))

# Regression
summary(Fit1 <- lm(Correct ~ Temperature, data = X))
summary(lm(abs(Fit1$residuals) ~ X$Temperature))

(tukey1 <- TukeyHSD(AOV1))
plot(tukey1)

# Create dummy codes
X$X1 <- X$Temperature == '50 deg.'
X$X2 <- X$Temperature == '90 deg.'

summary(zFit <- lm(Correct ~ X1 + X2, data=X))

# Deviate codes
X$d50 <- ifelse(X$Temperature == '50 deg.', 1, ifelse(X$Temperature == '90 deg.', -1, 0))
X$d70 <- ifelse(X$Temperature == '70 deg.', 1, ifelse(X$Temperature == '90 deg.', -1, 0))

summary(dFit <- lm(Correct~ d50 + d70, data=X))
TukeyHSD(aov(dFit))

# Orthogonal coding
tvHelmert <- matrix(c(-1/3, 2/3, 1/3, 1/2, 0, -1/2), ncol=2)

(contrasts(X$Temperature) <- tvHelmert)


############## Tues Nov 11th ##################

rm( list=ls() )

dfData <- data.frame( matrix( c(
1, 85, 100,
1, 80, 98,
1, 92, 105,
2, 86, 92,
2, 82, 99,
2, 95, 108,
3, 90, 95,
3, 87, 80,
3, 78, 82
), ncol=3, byrow=T ) )

colnames( dfData ) <- c( "zGroup", "PreTest", "PostTest" )

dfData$zGroup <- factor( dfData$zGroup )

dfData

# Model1 - this is simultaneous (BAD MODEL) -- order matters!
# This regressed our group out first
Model1 <- lm( PostTest ~ zGroup + PreTest, data=dfData )
# Inclu
summary( Model1 )
summary( aov( Model1 ) )

TukeyHSD( aov( Model1 ), 'zGroup' )
plot( TukeyHSD( aov( Model1 ), 'zGroup' ), las=1 )

# Model2 - hierarchical (GOOD MODEL) -- order matters!
# This regressed pretest out first
Model2 <- lm( PostTest ~ PreTest + zGroup, data=dfData )
summary( Model2 )
summary( aov( Model2 ) )

TukeyHSD( aov( Model2 ), 'zGroup' )
plot( TukeyHSD( aov( Model2 ), 'zGroup' ), las=1 )

############## Tues Dec 2nd ##################

rm( list=ls() )

load( '~/Documents/Classwork/Statistics_711/BtheB.RData' )

layout( matrix(1:2, nrow = 1 ) )

ylim <- range(BtheB[,grep("bdi", names(BtheB))], na.rm = TRUE)

boxplot(subset(BtheB, treatment == "TAU")[,grep("bdi", names(BtheB))],
        main = "Standard Treatment", ylab = "BDI", las=1,
        xlab = "Time (in months)", names = c(0, 2, 3, 5, 8), ylim = ylim)   

boxplot(subset(BtheB, treatment == "BtheB")[,grep("bdi", names(BtheB))],
        main = "Beat the Blues", ylab = "BDI", xlab = "Time (in months)",
        last=1, names = c(0, 2, 3, 5, 8), ylim = ylim)

BtheB$subject <- 1:nrow( BtheB )

dfLong <- reshape( data=BtheB, dir='long', idvar='subject', 
                   varying = c("bdi.2m", "bdi.3m", "bdi.5m", "bdi.8m"), times=c( 2, 3, 5, 8 ) )

dfLong <- dfLong[ order( dfLong$subject, dfLong$time ),
                  c( 'subject', 'drug', 'length', 'treatment', 'bdi.pre', 'time', 'bdi' )]

library( nlme )
library( multcomp )

# Random intercept only model
Fit1 <- lme( fixed = bdi ~ bdi.pre + time + treatment + drug, random = ~ 1 | subject,
   data = dfLong, na.action = na.omit, method='ML' )
   
# Random slope and intercept model (control for pre test) -- fixed effects
Fit2 <- lme( fixed = bdi ~ bdi.pre + time + treatment + drug + length, random = ~ 1 + time | subject,
   data = dfLong, na.action = na.omit, method='ML' )

# Compare the fit of the two models
anova( Fit1, Fit2 )

VarCorr(Fit1)
summary( Fit1 )

# ICC
( VarEsts <- as.numeric( VarCorr( Fit1 )[1:2 ] ) ) # vector of variance estimates
( ICC <- VarEsts[ 1 ] / sum( VarEsts ) ) # computes ICC

cmfit( Fit1 )

