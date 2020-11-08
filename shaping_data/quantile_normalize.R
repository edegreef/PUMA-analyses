#quick script to quantile-normalize data

library(BiocManager)
BiocManager::install("preprocessCore")
library(preprocessCore)

#load data file
data <- read.csv("firstbloom-greenness_87_chunks.csv", header=T)

#select columns that you want to quantile-normalize
sub <- data[,c(7,8)]

#remove NA's if necessary
sub <- na.omit(sub)

#transform into matrix
mat <- as.matrix(sub)

#quantile_normalize
qn <- normalize.quantiles(mat)

#prepping output to save as csv
library(tidyverse)
qn_named <- as.data.frame(qn)
qn_named <- select(qn_named, MODIS_old=V1, MODIS_new=V2)
write.csv(qn_named,"greenup_QN_87_chunks.csv")
