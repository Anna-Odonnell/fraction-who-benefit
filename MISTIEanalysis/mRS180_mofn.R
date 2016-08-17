##8/3/16

##Get the 95% CI for the lower bound (no restrictions, no baseline variable) 
##and the 95% CI for the upper bound (no restrictions, no baseline variable)
##using m-out-of-n bootstrap, when the outcome is 180-day mRS

library(boot)
library(lpSolveAPI)

rm(list=ls())
##########################################################################################################
##Functions##
##########################################################################################################
boundsNoCov_res <- function(ordinalScale, YT, YC, maxBen, maxHarm){
  nT <- length(YT) #number of treatment subjects
  nC <- length(YC) #number of control subjects
  
  if(nT == 0 | nC == 0){
    #stop("WARNING: YT or YC is empty")
    return(c(NA,NA,NA))
  } else {
    ordinalScale <- sort(ordinalScale, decreasing = FALSE)
    L <- length(ordinalScale)
    varCount <- L^2 #number of pi_{i,j}'s
    
    #The matrix of pi_{i,j}'s is a L x L matrix with varCount pi_{i,j}'s.
    
    scale <- nT * nC ##this can be made large to help with solving the linear program
    
    #############################################################################
    ##Calculate marginal cdf's##
    #############################################################################
    cdf_C <- rep(0,L)
    for (i in 1:L){
      cdf_C[i] <- sum(YC <= ordinalScale[i])/nC
    }
    cdf_T <- rep(0,L)
    for (i in 1:L){
      cdf_T[i] <- sum(YT <= ordinalScale[i])/nT
    }
    
    #############################################################################
    ##Which pi_{i,j}'s are affected by the restrictions?##
    #############################################################################
    restrictions <- matrix(0,nrow=L,ncol=L)
    for (i in 1:L){
      for (j in 1:L){
        if (ordinalScale[j]-ordinalScale[i]>maxBen | ordinalScale[i]-ordinalScale[j]>maxHarm){
          restrictions[i,j] <- 1
        }
      }
    }
    restrictions <- c(restrictions)
    #############################################################################
    ##LP for epsilon##
    #############################################################################
    lprec <- make.lp(0,(varCount+1))
    
    ##Setting objective function
    objfn <- c(rep(0,varCount),1)
    set.objfn(lprec,objfn)
    
    ##Setting non-negativity bounds
    set.bounds(lprec, lower = rep(0, (varCount+1)), upper = NULL)
    
    ##pi_{i,j}'s sum to 1
    add.constraint(lprec, xt = c(rep(1,varCount),0), "=", rhs = scale)
    
    ##incorporating the restrictions
    add.constraint(lprec, xt = c(restrictions,0), "=", rhs = 0)
    
    ##marginal cdf constraints
    pij.matrix <- matrix(1:varCount, nrow = L, ncol = L)
    
    for (i in 1:(L-1)){
      add.constraint(lprec, xt = c(rep(1,L*i),-1), "<=", rhs = (cdf_C[i]*scale), indices = c(pij.matrix[1:i,],(varCount+1))) 
      add.constraint(lprec, xt = rep(1,L*i+1), ">=", rhs = (cdf_C[i]*scale), indices = c(pij.matrix[1:i,],(varCount+1)))
    }
    
    for (i in 1:(L-1)){
      add.constraint(lprec, xt = c(rep(1,L*i),-1), "<=", rhs = (cdf_T[i]*scale), indices = c(pij.matrix[,1:i],(varCount+1))) 
      add.constraint(lprec, xt = rep(1,L*i+1), ">=", rhs = (cdf_T[i]*scale), indices = c(pij.matrix[,1:i],(varCount+1)))
    }
    
    ##Solving linear program
    eps.flag <- solve(lprec)
    
    if(eps.flag != 0){
      stop("WARNING: problem with LP for eps")
    }
    
    eps <- get.objective(lprec)/scale
    
    if (eps < 10 ^ (-10)){
      eps <- 0
    }
    
    rm(lprec)
    
    #############################################################################
    ##LP for lb and ub##
    #############################################################################
    
    lprec <- make.lp(0,varCount)
    
    ##Setting objective function
    objfn <- matrix(0, L, L)
    for (i in 1:L){
      for (j in 1:L){
        if (ordinalScale[j]>ordinalScale[i]){ 
          objfn[i,j] <- 1
        }
      }
    }
    objfn <- c(objfn)
    set.objfn(lprec,objfn)
    
    ##Setting non-negativity bounds 
    set.bounds(lprec, lower = rep(0, varCount), upper = NULL) 
    
    ##pi_{i,j}'s sum to 1
    add.constraint(lprec, xt = rep(1,varCount), "=", rhs = scale)
    
    ##incorporating the restrictions
    add.constraint(lprec, xt = restrictions, "=", rhs = 0)
    
    ##marginal cdf constraints
    for (i in 1:(L-1)){
      add.constraint(lprec, xt = rep(1,L*i), "<=", rhs = (cdf_C[i]+eps)*scale, indices = c(pij.matrix[1:i,])) 
      add.constraint(lprec, xt = rep(1,L*i), ">=", rhs = (cdf_C[i]-eps)*scale, indices = c(pij.matrix[1:i,]))
    }
    
    for (i in 1:(L-1)){
      add.constraint(lprec, xt = rep(1,L*i), "<=", rhs = (cdf_T[i]+eps)*scale, indices = c(pij.matrix[,1:i])) 
      add.constraint(lprec, xt = rep(1,L*i), ">=", rhs = (cdf_T[i]-eps)*scale, indices = c(pij.matrix[,1:i]))
    }
    
    lb.flag <- solve(lprec)
    
    if(lb.flag != 0){
      stop("WARNING: problem with LP for lb")
    }
    
    lb <- get.objective(lprec)/scale
    
    lp.control(lprec,sense='max')
    
    ub.flag <- solve(lprec)
    
    if(ub.flag != 0){
      stop("WARNING: problem with LP for ub")
    }
    
    ub <- get.objective(lprec)/scale
    
    return(c(lb,ub,eps))
  }
}


compute_estimators_on_bootstrap_replicate_of_dataset <- function(data,index,maxBen,maxHarm,subsampleSize){
  replicate_data_set <- data[sample(1:nrow(data),size=subsampleSize,replace=TRUE),]
  YC <- replicate_data_set[replicate_data_set$trt==0,]$y
  YT <- replicate_data_set[replicate_data_set$trt==1,]$y
  res <- boundsNoCov_res(ordinalScale, YT, YC, maxBen, maxHarm)
  return(res[1:2])
}

#############################################################################
##Get the dataset with 180-day mRS##
#############################################################################

setwd("~/Dropbox/research/data/MISTIEII")
load("MISTIEIIdata.Rdata")
##Each of these datasets (i.e., mRSData, clotData, baselineData) has 96 people (42 medical, 54 surgical)
##All inclusion criteria have been met except non-missing 180-day mRS score and baseline NIHSS score

##We only need mRSData and baselineData right now
rm(clotData)

##Also, we want to focus on the data at 180-days
mrs180 <- subset(mrsData, Follow_up_Visit == 180) 
rm(mrsData)
length(unique(mrs180$patientName)) ##96

##Merge mrs180 with baselineData
data <- merge(baselineData, mrs180, by = "patientName")

##Check that merge was successful
all(data$Group_Assigned.x == data$Group_Assigned.y) ##true

##Drop redundant columns
drops <- c("RunIn_Randomized.y","Group_Assigned.y")
data <- data[,!(names(data) %in% drops)]

rm(drops, mrs180, baselineData)

##Only include people with baseline NIHSS and 180-day mRS score
data$patientName[is.na(data$Enrollment_NIHSS_Total)] 
data <- subset(data, !is.na(Enrollment_NIHSS_Total))
nrow(data) ##95

sum(is.na(data$rankin_score)) ##6 people missing rankin_score
data$patientName[is.na(data$rankin_score)] 
data <- subset(data, !is.na(rankin_score)) 
nrow(data) ##89


control <- subset(data, Group_Assigned.x == "Medical") ##37
treatment <- subset(data, Group_Assigned.x == "Surgical") ##52

YC <- control$rankin_score 
YT <- treatment$rankin_score


##The scores are on a scale from 0 (no symptoms) to 6 (death)
##Recode the scores to be Levels 1 - 7 (higher levels being better):
##Level 1 = 6; Level 2 = 5; Level 3 = 4; Level 4 = 3; Level 5 = 2; Level 6 = 1; Level 7 = 0
levelfun <- function(x){ 
  if(x == 0){
    y = 7
  } else if (x == 1){
    y = 6
  } else if (x == 2){
    y = 5
  } else if (x == 3){
    y = 4
  } else if (x == 4){
    y = 3
  } else if (x == 5){
    y = 2
  } else {
    y = 1
  }
  return(y)
}

YC <- sapply(YC, levelfun)
YT <- sapply(YT, levelfun)

datasamp <- data.frame(trt = c(rep(1, length(YT)),rep(0, length(YC))), y = c(YT, YC))
rm(YC,YT, control, treatment, data)

#############################################################################
##Inputs##
#############################################################################
ordinalScale <- 1:7
maxBen <- 100
maxHarm <- 100
bootrep <- 5000
N <- nrow(datasamp)
q <- 0.95
k <- 500 

#############################################################################
##Get the m-out-of-n CI's for the lower and upper bounds##
#############################################################################

set.seed(252905)
m.candidates <- unique(ceiling(N * q^(0:k)))
m.candidates <- m.candidates[m.candidates>=20]
numCandidates <- length(m.candidates)

matLB <- matrix(NA, nrow = numCandidates, ncol = bootrep) 
matUB <- matrix(NA, nrow = numCandidates, ncol = bootrep)

for (l in 1:numCandidates){
  boot.obj <- boot(data = datasamp, statistic= compute_estimators_on_bootstrap_replicate_of_dataset, R=bootrep, maxBen = maxBen, maxHarm = maxHarm, subsampleSize = m.candidates[l])
  matLB[l,] <- t(sort(boot.obj$t[,1], na.last = TRUE))
  matUB[l,] <- t(sort(boot.obj$t[,2], na.last = TRUE))
}

matLB <- cbind(m.candidates, matLB)
matLB <- matLB[complete.cases(matLB),] ##get rid of candidate m's with any NA's (m was too small, nT or nC = 0) 
m.candidatesLB <- matLB[,1]
matLB <- matLB[,-1]

matUB <- cbind(m.candidates, matUB)
matUB <- matUB[complete.cases(matUB),] 
m.candidatesUB <- matUB[,1]
matUB <- matUB[,-1]

temp1 <- matLB[1:(nrow(matLB)-1),]-matLB[2:nrow(matLB),]
temp1 <- abs(temp1)
temp1 <- apply(temp1, 1, max)
m.chosenLB <- m.candidatesLB[min(which(temp1==min(temp1)))] 


temp1 <- matUB[1:(nrow(matUB)-1),]-matUB[2:nrow(matUB),]
temp1 <- abs(temp1)
temp1 <- apply(temp1, 1, max)
m.chosenUB <- m.candidatesUB[min(which(temp1==min(temp1)))]

boot.obj <- boot(data = datasamp, statistic= compute_estimators_on_bootstrap_replicate_of_dataset, R=10000, maxBen = maxBen, maxHarm = maxHarm, subsampleSize = m.chosenLB)
LB.CI_mn <- c(quantile(boot.obj$t[,1], probs = c(.025,0.975)), m.chosenLB)
##2.5%      97.5%            
##  0.0000000  0.3098081 81.0000000 

boot.obj <- boot(data = datasamp, statistic= compute_estimators_on_bootstrap_replicate_of_dataset, R=10000, maxBen = maxBen, maxHarm = maxHarm, subsampleSize = m.chosenUB)
UB.CI_mn <- c(quantile(boot.obj$t[,2], probs = c(.025,0.975)), m.chosenUB)
##2.5%      97.5%            
##  0.5000000  0.8518519 51.0000000 