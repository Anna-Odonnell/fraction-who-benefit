##5/20/16

##Plug-in Estimator for Binary (no harm) case

##Compute the bias and SE for the plug-in estimator, when it is available

###########################################################################################################
##Load bound estimates for n = 100, 500, 1000##
###########################################################################################################
setwd("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/allResults/binary_noHarm")

load("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/allResults/binary_noHarm/n100_mofn.Rdata")
bounds.n100 <- bounds

load("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/allResults/binary_noHarm/n500_mofn.Rdata")
bounds.n500 <- bounds

load("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/allResults/binary_noHarm/n1000_mofn.Rdata")
bounds.n1000 <- bounds

rm(bounds, LB.CI_mn, UB.CI_mn)


###########################################################################################################
##Function for computing standard deviation##
###########################################################################################################
SE.fun <- function(x){sqrt(mean((x-mean(x))^2))}

###########################################################################################################
##Proportion of simulations where plug-in estimator had a solution##
##Bias and SE, if plug-in estimator has a solution##
###########################################################################################################
bounds.n100 <- data.frame(bounds.n100)
names(bounds.n100) <- c("lower","upper","eps")
bounds.n500 <- data.frame(bounds.n500)
names(bounds.n500) <- c("lower","upper","eps")
bounds.n1000 <- data.frame(bounds.n1000)
names(bounds.n1000) <- c("lower","upper","eps")

##Proportion of simulations where plug-in estimator had a solution
mean(bounds.n100$eps>0)
mean(bounds.n500$eps>0)
mean(bounds.n1000$eps>0)


##Bias conditional on plug-in estimator having a solution
bounds.n100 <- subset(bounds.n100, eps == 0)
bounds.n500 <- subset(bounds.n500, eps == 0)
bounds.n1000 <- subset(bounds.n1000, eps == 0)

trueLB <- 0
trueUB <- 0

mean(bounds.n100$lower)-trueLB
mean(bounds.n100$upper)-trueUB
mean(bounds.n500$lower)-trueLB
mean(bounds.n500$upper)-trueUB
mean(bounds.n1000$lower)-trueLB
mean(bounds.n1000$upper)-trueUB

##Standard error conditional on plug-in estimator having a solution
SE.fun(bounds.n100$lower)
SE.fun(bounds.n100$upper)
SE.fun(bounds.n500$lower)
SE.fun(bounds.n500$upper)
SE.fun(bounds.n1000$lower)
SE.fun(bounds.n1000$upper)


