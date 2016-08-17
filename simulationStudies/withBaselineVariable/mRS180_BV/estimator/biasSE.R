##6/2/16

##Bias, standard error, and proportion of simulations with no bound estimates (i.e., there were no
##control or treatment subjects for some subgroup)

rm(list=ls())

##Load the true bounds
setwd("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/withBaselineVariable/mRS180_BV/boundParameters")
load("trueBounds.Rdata")

##These data frames will be used to store estimator properties, including bias and standard error
estProp2 <- data.frame(n = c(100,500,1000), biasLower = NA, biasUpper = NA, SElower = NA, SEupper = NA, propNA = 0)
estProp3 <- data.frame(n = c(100,500,1000), biasLower = NA, biasUpper = NA, SElower = NA, SEupper = NA, propNA = 0)
estProp4 <- data.frame(n = c(100,500,1000), biasLower = NA, biasUpper = NA, SElower = NA, SEupper = NA, propNA = 0)

##Function for getting SE
SE.fun <- function(x){sqrt(mean((x-mean(x))^2))}

##Function for getting estimator properties
estProp <- function(boundEst, trueBounds){
  propNA <- mean(is.na(bounds$eps))*100
  bounds.noNA <- subset(bounds, !is.na(eps))
  bias <- colMeans(bounds.noNA)-trueBounds
  SE <- c(SE.fun(bounds$lb), SE.fun(bounds$ub))
  return(c(bias[1:2], SE, propNA))
}

##Get estimator properties

setwd("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/withBaselineVariable/mRS180_BV/estimator")

load("subpop2_N100.Rdata")
estProp2[1,2:6] <- estProp(bounds, res2)
estProp2

load("subpop2_N500.Rdata")
estProp2[2,2:6] <- estProp(bounds, res2)
estProp2

load("subpop2_N1000.Rdata")
estProp2[3,2:6] <- estProp(bounds, res2)
estProp2


load("subpop3_N100.Rdata")
estProp3[1,2:6] <- estProp(bounds, res3)
estProp3

load("subpop3_N500.Rdata")
estProp3[2,2:6] <- estProp(bounds, res3)
estProp3

load("subpop3_N1000.Rdata")
estProp3[3,2:6] <- estProp(bounds, res3)
estProp3



load("subpop4_N100.Rdata")
estProp4[1,2:6] <- estProp(bounds, res4)
estProp4

load("subpop4_N500.Rdata")
estProp4[2,2:6] <- estProp(bounds, res4)
estProp4

load("subpop4_N1000.Rdata")
estProp4[3,2:6] <- estProp(bounds, res4)
estProp4



##Format the results into LaTeX tables
library(xtable)
table <- xtable(rbind(estProp2, estProp3, estProp4), align = "ccccccc")
digits(table) <- 3
print(table)
