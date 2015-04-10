# HOMEWORK 2
# Do partial and semipartial correlations of brain data

rm(list = ls())
library(ggplot2)
library(ppcor)

# DIRECTORY STUFF (tab delimited w. header)
dir = '~/Documents/Classwork/Statistics_711/HW1_2/brain.dat'
dat = read.table(dir, header = 1, sep = '\t')

r = cor.test(dat$FSIQ,dat$MRI_Count)

## Partial correlation partialling out gender:
pr = pcor.test(dat$FSIQ,dat$MRI_Count,dat$Female)
spr = spcor.test(dat$FSIQ,dat$MRI_Count,dat$Female)

## Do it manually:
# IQ
IQ.diff = dat$FSIQ - mean(dat$FSIQ)
IQ.var = IQ.diff^2
IQ.SS = sum(IQ.var)

# Brain Size
BS.diff = dat$MRI_Count - mean(dat$MRI_Count)
BS.var = BS.diff^2
BS.SS = sum(BS.var)

# Gender
G.diff = dat$Female - mean(dat$Female)
G.var = G.diff^2
G.SS = sum(G.var)

# Sum of products and correlation of IQ and brain size
SP.IQ_BS = sum(BS.diff * IQ.diff)
r.IQ_BS = SP.IQ_BS / sqrt(IQ.SS*BS.SS)

# Sum of products and correlation of IQ and gender
SP.IQ_G = sum(IQ.diff * G.diff)
r.IQ_G = SP.IQ_G / sqrt(IQ.SS * G.SS)

# Sum of products and correlation of brain size and gender
SP.BS_G = sum(BS.diff * G.diff)
r.BS_G = SP.BS_G / sqrt(BS.SS * G.SS)

## Partial correlation, holding gender constant:
(r.IQ_BS - (r.IQ_G*r.BS_G))/(sqrt((1-r.IQ_G^2) * (1-r.BS_G^2)))
## Semi-partial correlation, holding gender constant:
(r.IQ_BS - (r.IQ_G*r.BS_G))/sqrt((1-r.BS_G^2))








