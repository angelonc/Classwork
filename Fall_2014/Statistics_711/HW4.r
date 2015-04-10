rm(list = ls())

x <- data.frame(matrix(c(
    1, 1, 24,
    1, 1, 26,
    1, 1, 25,
    1, 1, 24,
    1, 1, 27,
    1, 1, 24,
    1, 1, 27,
    1, 1, 23,
    1, 0, 15,
    1, 0, 17,
    1, 0, 20,
    1, 0, 16,
    0, 1, 25,
    0, 1, 29,
    0, 1, 27,
    0, 0, 19,
    0, 0, 18,
    0, 0, 21,
    0, 0, 20,
    0, 0, 21,
    0, 0, 22,
    0, 0, 19 ), ncol = 3, byrow = T))

# Factors etc.
names(x)    <- c('Sex', 'Education', 'Income')
x$Sex       <- factor(x$Sex, labels = c('Male','Female'))
x$Education <- factor(x$Education, labels = c('No College', 'College'))

# Table sums
table(Sex = x$Sex, Education = x$Education)

# Display means and standard deviations
tapply(x$Income, list(x$Sex, x$Education), mean)
tapply(x$Income, list(x$Sex, x$Education), sd)

# Tells R to make contrasts sum to 0 by default
options(contrasts = c(unordered = 'contr.sum', ordered = 'contr.poly')) # What does this do?

# Initial model
(aFit1 <- aov(Income ~ Sex * Education, data = x))
summary(aFit1)

# Switch the variables - this does have an impact... gender has a significant main effect
# These use type I, which we don't care about
(aFit2 <- aov(Income ~ Education * Sex, data = x))
summary(aFit2)

# Multiple comparisons
TukeyHSD(aFit2, which = 'Education')
TukeyHSD(aFit2, which = 'Sex')


# Type III Sums of Squares (simultaneous effects) - gender comes out significant again
# Good for balanced designs
drop1(aFit1, .~., test = 'F')

# Car library
library(car)

Anova(aFit1, type = 2)
Anova(aFit2, type = 3)

library(ggplot2)
library(plyr)

Means <- ddply(x, .(Education, Sex), summarise, val = mean(Income))

ggplot(x, aes(x = Education, y = Income, color = Sex)) +
    geom_point(data = Means, aes(y = val)) +
    geom_boxplot() +
    geom_line(data = Means, aes(y = val, group = Sex)) +
    theme_bw()

# Effect sizes
# R-squared (overestimated relationships with small samples)
# Omega-squared
# Cohen's d


# Additional hypothesis - regressing out educational level, is there an effect of gender on income
GenderOnly <- lm(Income ~ Education, data = x)
summary(lm(Income ~ GenderOnly$resid, data = x))
#Anova(GenderOnly$resid, type = 2)
