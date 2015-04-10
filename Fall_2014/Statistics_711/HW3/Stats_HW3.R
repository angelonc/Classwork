# HOMEWORK 3 TASKS:
# Analyze regression data
# WRite up results
# Think about whether data satisfy assumptions of linear model
# Be prepared to talk about our healthcare related behaviors
# Consequences of what you observe in your data when diagnosing your model

rm(list = ls())

library(ggplot2)  # Plotting
library(car)     # Diagnosing
library(GGally)  # Nice plots

dat = read.table("~/Documents/Classwork/Statistics_711/HW3/regress.dat", header=1)
head(dat)

# DATA FIELDS:
# timedrs = number of health care visits
# phyheal = number of physical health symptoms
# menheal = number of mental health symptoms
# stress  = stress from acute life changes

tvDescriptive <- function( x ){
    z <- as.numeric( x )
    N <- sum( !is.na( x ) )
    NMiss <- sum( is.na( x ) )
    Mean <- mean( z )
    Median <- median( z )
    P25 <- quantile( z, .25 )
    P75 <- quantile( z, .75 )
    SD <- sd( z )
    Min <- min( z )
    Max <- max( z )
    SE <- SD / N^.5
    CV <- SD / Mean
    LCL <- Mean - SE * qt( .975, N - 1 )
    UCL <- Mean + SE * qt( .975, N - 1 )
    
    M3 <- sum( ( z - Mean )^3 ) # This is the 3rd moment about the mean
    Skew <- M3 / ( ( N - 1 ) * SD^3 )
    
    M4 <- sum( ( z - Mean )^4 ) # This is the 4th moment about the mean
    Kurtosis <- M4 / ( ( N - 1 ) * SD^4 )
    
    return( c( N=N, Missing=NMiss, Mean=Mean, Median=Median, 'P25'=P25, 'P75'=P75,
              SD=SD, SE=SE, Min=Min, Max=Max, CV=CV, Skew=Skew,
              Kurtosis=Kurtosis, LCL=LCL, UCL=UCL ) )
}

# Apply the function we defined above to our data then transposing
t(apply(dat[,2:5],2,tvDescriptive))

cor(dat[,2:5])

ggpairs(dat[,2:5])

# DATA DIAGNOSES:
# VIF - this uses matrix algebra
# Solve: find inverse of correlation matrix
# Diag: returns diagonal of a square matrix as a vector
diag( solve( cor( dat[ , c( 'phyheal', 'menheal', 'stress' ) ] ) ) )

# Linear model:
summary(Fit1 <- lm(timedrs ~ phyheal + menheal + stress, data=dat))

confint(Fit1)

vif(Fit1)

plotDat <- data.frame(timedrs=dat$timedrs, fitted.values = Fit1$fitted.values, residuals = Fit1$residuals)

# Observed vs Expected scatter plot
ggplot(plotDat, aes( x=fitted.values, y=timedrs ) ) +
geom_point(shape=19) + # Use filled circles
geom_smooth(method=lm, # Add linear regression line
se=FALSE) # Don't add shaded confidence region

# Residual Plot
ggplot(plotDat, aes( x=fitted.values, y=residuals ) ) +
geom_point(shape=19) + # Use filled circles
geom_smooth(method=lm, # Add linear regression line
se=FALSE) # Don't add shaded confidence region

# Standardize the residuals
sResiduals <- rstandard( Fit1 )
hist( sResiduals, freq=F, las=1, xlab='Standardized Residuals', main='Histogram of Residuals', col='navy'
)
curve( dnorm, add=T, col='red', lwd=3 )
plot( Fit1 )

# Identify observations that have D > 4 / (n - k - 1)
cutoff = 4 / (nrow(dat) - 3 - 1)
plot(Fit1, which=4, cook.levels = cutoff, las = 1)
abline(h=cutoff)

# Z-scored all data
d <- apply( dat[ , 2:5 ], 2, FUN=function( x ) x - mean( x ) )
summary( lm( timedrs ~ phyheal + menheal + stress, data=data.frame( d ) ) )

library(MASS)

boxcox(Fit1)



graphics.off()
