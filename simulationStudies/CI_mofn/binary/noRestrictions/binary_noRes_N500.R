rm(list=ls())

temp <- commandArgs(TRUE)
seed <- as.numeric(temp[1])
 

#############################################################################
##Libraries##
#############################################################################
library(lpSolveAPI)
library(boot)
#############################################################################
##Function for estimating bounds##
#############################################################################
##This function calculates the bound estimates.##
##Assumes an ordinal outcome. The higher the ranking, the better.
##YT and YC are the vector of outcomes for the treatment and control groups.
##maxBen is the maximum possible benefit (i.e., YT - YC <= maxBen)
##maxHarm is the maximum possible harm (i.e., YC - YT <= maxHarm)

boundsNoCov_res <- function(YT, YC, maxBen, maxHarm){
  nT <- length(YT) #number of treatment subjects
  nC <- length(YC) #number of control subjects
  
  if(nT == 0 | nC == 0){
    return(c(NA,NA,NA))
  } else {
    YT_sort <- sort(unique(YT))
    YC_sort <- sort(unique(YC))
    
    mT <- length(YT_sort)
    mC <- length(YC_sort)
    
    if (mT == 1 & mC == 1){
      if ((YT_sort[1]-YC_sort[1]>maxBen) | (YC_sort[1]-YT_sort[1]>maxHarm)){
        return(c(NA,NA,NA)) ##the restrictions prevent the sum of the probabilities from being 1
      } else {
        temp1 <- as.numeric(YT_sort[1] > YC_sort[1]) 
        return(c(temp1,temp1,0))
      }
    } else {
      varCount <- mT * mC #number of pi_{i,j}'s
      
      #The matrix of pi_{i,j}'s is a mC x mT matrix with varCount pi_{i,j}'s.
      #The rows of the matrix are labelled by YC_sort, and the columns by YT_sort.
      
      scale <- max(nT,nC)
      
      #############################################################################
      ##Calculate marginal cdf's##
      #############################################################################
      cdf_C <- rep(0,mC)
      for (i in 1:mC){
        cdf_C[i] <- sum(YC <= YC_sort[i])/nC
      }
      cdf_T <- rep(0,mT)
      for (i in 1:mT){
        cdf_T[i] <- sum(YT <= YT_sort[i])/nT
      }
      
      #############################################################################
      ##Which pi_{i,j}'s are affected by the restrictions?##
      #############################################################################
      restrictions <- matrix(0,nrow=mC,ncol=mT)
      for (i in 1:mC){
        for (j in 1:mT){
          if (YT_sort[j]-YC_sort[i]>maxBen | YC_sort[i]-YT_sort[j]>maxHarm){
            restrictions[i,j] <- 1
          }
        }
      }
      restrictions <- c(restrictions)
      #############################################################################
      ##LP for epsilon##
      #############################################################################
      lprec <- make.lp(0,(varCount+1))
      
      #print("setting objective function: 1")
      
      ##Setting objective function
      objfn <- c(rep(0,varCount),1)
      set.objfn(lprec,objfn)
      
      #print("setting constraints: 1")
      
      ##Setting non-negativity bounds
      set.bounds(lprec, lower = rep(0, (varCount+1)), upper = NULL)
      
      ##pi_{i,j}'s sum to 1
      add.constraint(lprec, xt = c(rep(1,varCount),0), "=", rhs = scale)
      
      ##incorporating the restrictions
      add.constraint(lprec, xt = c(restrictions,0), "=", rhs = 0)
      
      ##marginal cdf constraints
      pij.matrix <- matrix(1:varCount, nrow = mC, ncol = mT)
      
      if(mC>1){
        for (i in 1:(mC-1)){
          add.constraint(lprec, xt = c(rep(1,mT*i),-1), "<=", rhs = (cdf_C[i]*scale), indices = c(pij.matrix[1:i,],(varCount+1))) 
          add.constraint(lprec, xt = rep(1,mT*i+1), ">=", rhs = (cdf_C[i]*scale), indices = c(pij.matrix[1:i,],(varCount+1)))
        }
      } 
      
      if(mT>1){
        for (i in 1:(mT-1)){
          add.constraint(lprec, xt = c(rep(1,mC*i),-1), "<=", rhs = (cdf_T[i]*scale), indices = c(pij.matrix[,1:i],(varCount+1))) 
          add.constraint(lprec, xt = rep(1,mC*i+1), ">=", rhs = (cdf_T[i]*scale), indices = c(pij.matrix[,1:i],(varCount+1)))
        }
      } 
      
      #print("solving LP's: 1")
      
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
      
      #print("setting objective function: 2")
      
      ##Setting objective function
      objfn <- matrix(0, mC, mT)
      for (i in 1:mC){
        for (j in 1:mT){
          if (YT_sort[j]>YC_sort[i]){ 
            objfn[i,j] <- 1
          }
        }
      }
      objfn <- c(objfn)
      set.objfn(lprec,objfn)
      
      #print("setting constraints: 2")
      
      ##Setting non-negativity bounds 
      set.bounds(lprec, lower = rep(0, varCount), upper = NULL) 
      
      ##pi_{i,j}'s sum to 1
      add.constraint(lprec, xt = rep(1,varCount), "=", rhs = scale)
      
      ##incorporating the restrictions
      add.constraint(lprec, xt = restrictions, "=", rhs = 0)
      
      ##marginal cdf constraints
      if(mC>1){
        for (i in 1:(mC-1)){
          add.constraint(lprec, xt = rep(1,mT*i), "<=", rhs = (cdf_C[i]+eps)*scale, indices = c(pij.matrix[1:i,])) 
          add.constraint(lprec, xt = rep(1,mT*i), ">=", rhs = (cdf_C[i]-eps)*scale, indices = c(pij.matrix[1:i,]))
        }
      }
      
      if(mT>1){
        for (i in 1:(mT-1)){
          add.constraint(lprec, xt = rep(1,mC*i), "<=", rhs = (cdf_T[i]+eps)*scale, indices = c(pij.matrix[,1:i])) 
          add.constraint(lprec, xt = rep(1,mC*i), ">=", rhs = (cdf_T[i]-eps)*scale, indices = c(pij.matrix[,1:i]))
        }
      }
      
      #print("solving LP's: 2")
      
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
}



#############################################################################
##Choose m##
#############################################################################

compute_estimators_on_bootstrap_replicate_of_dataset <- function(data,index,maxBen,maxHarm,subsampleSize){
  replicate_data_set <- data[sample(1:nrow(data),size=subsampleSize,replace=TRUE),]
  YC <- replicate_data_set[replicate_data_set$trt==0,]$y
  YT <- replicate_data_set[replicate_data_set$trt==1,]$y
  res <- boundsNoCov_res(YT, YC, maxBen, maxHarm)
  return(res[1:2])
}

#############################################################################
##Design the true population##
#############################################################################
YCpop <- c(0,1)
YTpop <- c(0,1)
maxBen <- 100
maxHarm <- 100

#############################################################################
##Simulate randomized trials##
#############################################################################
nsim <- 1
bootrep <- 5000
N <- 500 ##update
q <- 0.95
k <- 500 ##update

m.candidates <- unique(ceiling(N * q^(0:k)))
m.candidates <- m.candidates[m.candidates>=20]
numCandidates <- length(m.candidates)

bounds <- matrix(data = NA, nrow = nsim, ncol = 3) #bounds and eps
LB.CI_mn <- matrix(data = NA, nrow = nsim, ncol = 3) ##first two elements are limits of CI, third element is the m
UB.CI_mn <- matrix(data = NA, nrow = nsim, ncol = 3)

for (s in 1:nsim){
  set.seed(seed)
  YCsamp <- sample(YCpop, size = N/2, replace = TRUE)
  YTsamp <- sample(YTpop, size = N/2, replace = TRUE)
  
  bounds[s,] <- boundsNoCov_res(YTsamp, YCsamp, maxBen, maxHarm)
  
  datasamp <- data.frame(trt=c(rep(1,length(YTsamp)),rep(0,length(YCsamp))),y=c(YTsamp,YCsamp)) #trt is 1 if treatment and 0 if control, y is the outcome
  
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
  temp1 <- boot.obj$t[,1]
  LB.CI_mn[s,] <- c(quantile(temp1, probs = c(.025,0.975)), m.chosenLB)
  
  boot.obj <- boot(data = datasamp, statistic= compute_estimators_on_bootstrap_replicate_of_dataset, R=10000, maxBen = maxBen, maxHarm = maxHarm, subsampleSize = m.chosenUB)
  temp1 <- boot.obj$t[,2]
  UB.CI_mn[s,] <- c(quantile(temp1, probs = c(.025,0.975)), m.chosenUB)
}

save(bounds, LB.CI_mn, UB.CI_mn, file=paste("binary500_noRes","-","seed", seed, "-", format(Sys.time(), "%Y%m%d-%H%M"),".Rdata", sep=""))
