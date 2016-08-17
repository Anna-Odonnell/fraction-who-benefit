##8/7/16

rm(list=ls())

##MAKE SURE TO DO THIS: Run the code in the file boundFunctions.R in the folder entitled "boundFunctions"##

setwd("~/Dropbox/research/github/fraction-who-benefit/demo")

##Load the simulated dataset
load("simulatedDataset.Rdata")


##The dataset is in wide format. There are 300 participants, one row per participant.
##A = treatment assignment (0 if control, 1 if treatment)
##X = categorical baseline variable
##Y = ordinal outcome, with larger values corresponding to better outcomes

##In this case, the ordinal outcome Y has four possible levels (1 - 4).
ordinalScale <- 1:4

YT <- data$Y[data$A == 1] ##outcomes of participants assigned to treatment
YC <- data$Y[data$A == 0] ##outcomes of participants assigned to control

YT1 <- data$Y[data$A == 1 & data$X == 1] ##outcomes of participants in the first stratum of X who are assigned to treatment 
YC1 <- data$Y[data$A == 0 & data$X == 1] ##outcomes of participants in the first stratum of X who are assigned to control
YT2 <- data$Y[data$A == 1 & data$X == 2] ##outcomes of participants in the second stratum of X who are assigned to treatment
YC2 <- data$Y[data$A == 0 & data$X == 2] ##outcomes of participants in the second stratum of X who are assigned to control

##The baseline variable X has two categories labeled 1 and 2.



######################################################################################################
##Compute bound estimates##
######################################################################################################

##EXAMPLES OF BOUND ESTIMATES, WITHOUT BASELINE VARIABLE##

##Notes:
##Each result is a vector with three elements. The first element is the lower bound estimate,
##the second is the upper bound estimate, and the third is the value of epsilon_bar.
##If you do not want to impose restrictions, set the inputs maxBen and maxVal to be a value
##larger than the span of the ordinal scale.

##Without restrictions, without baseline variable
est1 <- boundsNoCov_res(ordinalScale, YT, YC, maxBen = 100, maxHarm = 100)

##With no harm assumption, without baseline variable
est2 <- boundsNoCov_res(ordinalScale, YT, YC, maxBen = 100, maxHarm = 0)

##With assumption that benefit is at most 2 levels, without baseline variable
est3 <- boundsNoCov_res(ordinalScale, YT, YC, maxBen = 2, maxHarm = 100)

##With no harm assumption and the assumption that benefit is at most 1 level,without baseline variable
est4 <- boundsNoCov_res(ordinalScale, YT, YC, maxBen = 1, maxHarm = 0)



##EXAMPLES OF BOUND ESTIMATES, WITH BASELINE VARIABLE##

##Notes:
##Each result is a vector with eight elements. Elements 1 and 2
##are the lower and upper bound estimates for the population.
##Elements 3-5 are the lower bound estimate, the upper estimate, and value of epsilon_bar_1
##for the subpopulation with X = 1. 
##Elements 6-8 are the lower bound estimate, the upper estimate, and value of epsilon_bar_2
##for the subpopulation with X = 2.
##Our method can handle a baseline variable with any finite number of
##categories (it does not need to be binary). The code we provide 
##(in the boundFunctions folder on GitHub) and the examples below incorporate a binary X, 
##but it can be extended to handle X with more than two categories. 

##Without restrictions, with baseline variable
est5 <- boundsCov_res(ordinalScale, YT1, YC1, YT2, YC2, maxBen = 100, maxHarm = 100)

##With no harm assumption, with baseline variable
est6 <- boundsCov_res(ordinalScale, YT1, YC1, YT2, YC2, maxBen = 100, maxHarm = 0)

##With assumption that benefit is at most 2 levels, with baseline variable
est7 <- boundsCov_res(ordinalScale, YT1, YC1, YT2, YC2, maxBen = 2, maxHarm = 100)

##With no harm assumption and assumption that benefit is at most 1 level, with baseline variable
est8 <- boundsCov_res(ordinalScale, YT1, YC1, YT2, YC2, maxBen = 1, maxHarm = 0)

######################################################################################################
##Get two-sided 95% CI for lower bound and two-sided 95% CI for upper bound,
##using m-out-of-n bootstrap (without restrictions, without baseline variable)
######################################################################################################

##MAKE SURE TO DO THIS: Run the code in the file mOfn.R in the folder entitled "boundFunctions"

set.seed(35841358)

N <- nrow(data) ##sample size of dataset

maxBen <- 100 ##indicate that no restrictions are to be imposed
maxHarm <- 100

q <- 0.95 ##inputs for determining candidate values for m 
k <- 500 
m.candidates <- unique(ceiling(N * q^(0:k)))
m.candidates <- m.candidates[m.candidates>=20] ##these are the candidate values for m

bootrep <- 5000 ##This is the number of replicated datasets to be generated
                ##for each candidate value m. Let it be a large number.
bootrep.chosenM <- 10000 ##This is the number of replicated datasets to be generated
                         ##for the value m that is selected. 


##STEP 1: Selecting m for the lower bound and separately for the upper bound##
numCandidates <- length(m.candidates)
  
matLB <- matrix(NA, nrow = numCandidates, ncol = bootrep) 
matUB <- matrix(NA, nrow = numCandidates, ncol = bootrep)


for (l in 1:numCandidates){
  boot.obj <- boot(data = data, statistic= compute_estimators_on_bootstrap_replicate_of_dataset, R=bootrep, maxBen = maxBen, maxHarm = maxHarm, subsampleSize = m.candidates[l])
  matLB[l,] <- t(sort(boot.obj$t[,1], na.last = TRUE))
  matUB[l,] <- t(sort(boot.obj$t[,2], na.last = TRUE))
}
##This loop gets the lower and upper bound estimates of the 5,000 bootstrap replicated datasets. The
##results are stored in matLB and matUB, respectively. Since this is done for each candidate value m,
##it is the most time-consuming step of m-out-of-n bootstrap. On my laptop, lines 117 to 121 took
##12 minutes to run. It can be sped up through parallelization, since the results for the different
##candidates can be computed independently.

##For each candidate value, there is a row in matLB (matUB) with the lower (upper) bound estimates for that candidate value.
##The lower (upper) bound estimates are ordered from smallest to largest, with any NA's at the end.

##Now we select the value m for the lower bound
matLB <- cbind(m.candidates, matLB) ##Label each row with the candidate value to which it corresponds.
matLB <- matLB[complete.cases(matLB),] ##get rid of candidate m's with any NA's (m was too small, nT or nC = 0) 
m.candidatesLB <- matLB[,1] ##these are the candidates remaining after excluding those with NA's
matLB <- matLB[,-1] ##Take off the labels. We don't want these when choosing m.

temp1 <- matLB[1:(nrow(matLB)-1),]-matLB[2:nrow(matLB),] ##Take the difference between adjacent rows (row 1 - row 2, row2 - row3, row3 - row4, ...)
temp1 <- abs(temp1) ##Take the absolute value (|row1 - row2|, |row2 - row3|, |row3 - row4|, ...)
temp1 <- apply(temp1, 1, max) ##Find the maximum (max{|row1 - row2|}, max{|row2 - row3|}, max{|row3 - row4}, ...)

m.chosenLB <- m.candidatesLB[min(which(temp1==min(temp1)))] 
##Choose the value m with the smallest max absolute difference. 
##If multiple m attain the smallest max absolute difference, choose the largest one.
##m.chosenLB is the value m selected for lower bound

##Separately, but using the same procedure, we select the value m for the upper bound
matUB <- cbind(m.candidates, matUB)
matUB <- matUB[complete.cases(matUB),] 
m.candidatesUB <- matUB[,1]
matUB <- matUB[,-1]
  
temp1 <- matUB[1:(nrow(matUB)-1),]-matUB[2:nrow(matUB),]
temp1 <- abs(temp1)
temp1 <- apply(temp1, 1, max)

m.chosenUB <- m.candidatesUB[min(which(temp1==min(temp1)))] 
##m.chosenUB is the value m selected for upper bound

##STEP 2: Run bootstrap using selected values of m##
boot.obj <- boot(data = data, statistic= compute_estimators_on_bootstrap_replicate_of_dataset, R=bootrep.chosenM, maxBen = maxBen, maxHarm = maxHarm, subsampleSize = m.chosenLB)
LB.CI_mn <- c(quantile(boot.obj$t[,1], probs = c(.025,0.975)), m.chosenLB)
##On my laptop, lines 160-161 took 27 seconds to run.


boot.obj <- boot(data = data, statistic= compute_estimators_on_bootstrap_replicate_of_dataset, R=bootrep.chosenM , maxBen = maxBen, maxHarm = maxHarm, subsampleSize = m.chosenUB)
UB.CI_mn <- c(quantile(boot.obj$t[,2], probs = c(.025,0.975)), m.chosenUB)
##On my laptop, lines 165-166 took 27 seconds to run.

##STEP 3: Gather the results in a data frame and print the results##
result <- rbind(LB.CI_mn, UB.CI_mn)
result <- data.frame(result)
names(result) <- c("lowerLimit", "upperLimit", "m")
row.names(result) <- c("lower bound", "upper bound")

result ##This gives the CI's


######################################################################################################
##Review all results from the code in this R script##
######################################################################################################

est1
est2
est3
est4
est5
est6
est7
est8
result

##Refer to Supplementary Materials for the output