rm(list = ls())

# Make group mats (first column = IQ, 2nd = score)
perm = cbind(c(132,122,118,121,101,117,122,122,119,126,113,122,111,93,124,129,131,112,112,115,123),  c(14,10,15,10,14,11,5,8,22,25,25,22,18,17,20,17,29,17,20,21,17))
auth = cbind(c(138,115,112,127,122,120,111,119,126,110,115,116,127,116,122,123,125,114,129), c(28,13,13,13,10,18,13,13,13,19,14,20,17,17,13,21,16,27,34))
demo = cbind(c(94,114,120,102,110,105,119,120,127,112,128,110,114,120,118,111), c(9,18,17,13,23,21,25,24,22,21,24,20,21,24,25,19))

# Assign to data structure
dat = cbind(c(rep(1,dim(perm)[1]),rep(2,dim(auth)[1]),rep(3,dim(demo)[1])),rbind(perm,auth,demo))
dat = data.frame(dat)
colnames(dat) = c('tGroup','IQ','Score')
dat$tGroup = factor(dat$tGroup, labels = c('Perm','Auth','Demo'))

# Do ANCOVA to see if teaching style is associated with test performace keeping IQ constant
# First regress out IQ and include interaction
Fit1 <- lm(Score ~ IQ + tGroup + IQ:tGroup, data = dat)
# ANOVA
anova(Fit1)

# With only the two predictors.. no interaction so we are justified to run this
Fit2 <- lm(Score ~ IQ + tGroup, data = dat)
anova(Fit2)

# Extract regression coefficients intercepts
coeffs = coef(Fit2)
B0.Permissive = coeffs[1]
B2.Authoritarian = coeffs[2] + B0.Permissive
B3.Democratic = coeffs[3] + B0.Permissive
slopeCommon = coeffs[2]

# Plot different conditions
plot(perm, xlab = 'IQ', ylab = 'Score', pch = 19, col = 'green')
points(auth, pch = 19, col = 'red')
points(demo, pch = 19, col = 'blue')

legend(x = 'topleft', legend = levels(data$tGroup), pch = c(19,19,19), col = c('green','red','blue'))
abline(B0.Permissive, slopeCommon, col='green')
abline(B2.Authoritarian, slopeCommon, col='red')
abline(B3.Democratic, slopeCommon, col='blue')
