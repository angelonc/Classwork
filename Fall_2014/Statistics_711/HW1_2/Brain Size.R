# Brain Size and Mental Capacity
# T. Victor

# Clear memory
rm( list=ls() )

library( ggplot2 )

dir <- "~/Dropbox/Today/"

# Read the brain data from the file
dfBrain <- read.table( paste( dir, "brain.dat", sep="" ), header=TRUE )

# Generate a scatter plot of the data
ggplot(dfBrain, aes( x=MRI_Count, y=FSIQ ) ) +
   xlab( 'Brain Size (MRI Count in Voxels)' ) +
   ylab( 'Full Scale Intelligence Quotient (Wechsler)' ) +
   geom_point( ) + # Points for the scatter plot
   geom_smooth( method=lm, se=FALSE )   # Add linear regression line

# tell R that the female variable is a factor
dfBrain$Female <- factor( dfBrain$Female, labels=c("Male", "Female"))

( Mean.IQ <- mean( dfBrain$FSIQ ) )
( Mean.BS <- mean( dfBrain$MRI_Count ) )

dfBrain$IQ.d <- dfBrain$FSIQ - Mean.IQ
dfBrain$IQ.d.2 <- dfBrain$IQ.d^2

dfBrain$BS.d <- dfBrain$MRI_Count - Mean.BS
dfBrain$BS.d.2     <- dfBrain$BS.d^2

dfBrain$Prod <- dfBrain$IQ.d * dfBrain$BS.d

( SP <- sum( dfBrain$Prod ) )
( SS.IQ <- sum( dfBrain$IQ.d.2 ) )
( SS.BS <- sum( dfBrain$BS.d.2 ) )

SP / sqrt( SS.IQ * SS.BS )

# return descriptive statistics for the quantitative variables
apply( dfBrain[ , 2:7 ], 2, function( x ) c(
                                        Mean    = mean( x, na.rm=TRUE )
                                       ,Variance = var( x, na.rm=TRUE )
                                       ,SD      = sd( x, na.rm=TRUE )
                                       ,Minimum = min( x, na.rm=TRUE )
                                       ,Maximum = max( x, na.rm=TRUE )
                                       ,Median  = median( x, na.rm=TRUE )
                                       ,Valid   = sum( !is.na( x ) )
                                       ,Missing = sum( is.na( x ) )
                                      )
)

# the correlation between brain size and full-scale IQ
cor.test( dfBrain$FSIQ, dfBrain$MRI_Count )

# the correlation between brain size and full-scale IQ by sex
by( dfBrain, dfBrain$Female, function(x) cor.test( x$FSIQ, x$MRI_Count ) )

ggplot(dfBrain, aes( x=MRI_Count, y=FSIQ, color=Female ) ) +
   xlab( 'Brain Size (MRI Count in Voxels)' ) +
   ylab( 'Full Scale Intelligence Quotient (Wechsler)' ) +
   geom_point( ) + # Points for the scatter plot
   geom_smooth( method=lm, fillrange=TRUE, se=FALSE )
   
   
