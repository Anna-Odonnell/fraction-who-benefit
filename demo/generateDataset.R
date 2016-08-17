##8/7/16

########################################################################################################
##Simulate a trial##
########################################################################################################
rm(list=ls())

set.seed(567698290)

n <- 300 
data <- data.frame(A = sample(0:1, size = n, replace = TRUE), X = sample(1:2, size = n, replace = TRUE))
data$Y[data$A == 0] <- sample(1:4, size = sum(data$A == 0), replace = TRUE, prob = c(1/4, 1/4, 1/4, 1/4)) 
data$Y[data$A == 1 & data$X == 1] <- sample(1:4, size = sum(data$A == 1 & data$X == 1), replace = TRUE, prob = c(1/8, 1/8, 1/2, 1/4))  
data$Y[data$A == 1 & data$X == 2] <- sample(1:4, size = sum(data$A == 1 & data$X == 2), replace = TRUE, prob = c(1/4, 1/4, 1/4, 1/4))  

setwd("~/Dropbox/research/github/fraction-who-benefit/demo")
save(data, file = "simulatedDataset.Rdata")
