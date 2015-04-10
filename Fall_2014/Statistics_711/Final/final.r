rm(list = ls())

# Load data
dat = read.table("~/Documents/Classwork/Statistics_711/Final/final_exam_data.csv", header=1, sep=',')

# Factors/labels
names(dat) = c('Treatment', 'Setting', 'Well_Being')
dat$Treatment = factor(dat$Treatment, labels = c('Freud','CBT','Fam_Sys'))
dat$Setting = factor(dat$Setting, labels = c('Home','Group','One_One'))
dat = data.frame(dat)
head(dat)

# Summary stats
hist(dat$Well_Being, main = "Well-Being Frequency Distribution", xlab = "Well-Being")
group.means <- aggregate( Well_Being ~ Setting + Treatment, dat, mean)
interaction.plot(dat$Setting, dat$Treatment, dat$Well_Being, fun = mean, xlab = "Setting", trace.label = "Treatment", ylab = "Average Well-Being", main = "Group Means by Condition")

aggregate( Well_Being ~ Setting + Treatment, dat, sd)
aggregate( Well_Being ~ Setting + Treatment, dat, sum)

# Set ANOVA options for car - ensures that contrast sum to 0
options(contrasts = c(unordered = 'contr.sum', ordered = 'contr.poly'))

library(ggplot2)
library(car)

# Initial model parameters (Type 1 SS)
aov.fit <- aov(Well_Being ~ Treatment * Setting, data = dat)

# Use car package to get type 2
fit1 = Anova(aov.fit, type = 2)
eta.squared.int <- fit1[[1]][[3]]/sum(fit1[[1]])

testInteractions(aov.fit)
testInteractions(aov.fit, pairwise="Treatment", across="Setting")
testInteractions(aov.fit, pairwise="Setting", across="Treatment")


# ANOVA 2
aov.fit2 <- aov(dat$Well_Being[dat$Treatment=="Freud"] ~ dat$Setting[dat$Treatment=="Freud"])
fit2 <- Anova(aov.fit2, type = 2)
testInteractions(aov.fit2)
