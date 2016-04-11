rm(list=ls())

#############################################################################
##Libraries##
#############################################################################
library(lpSolveAPI)
library(boot)
#############################################################################
##Function for estimating bounds##
#############################################################################
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


#############################################################################
##Choose m##
#############################################################################

compute_estimators_on_bootstrap_replicate_of_dataset <- function(data,index,maxBen,maxHarm,subsampleSize){
  replicate_data_set <- data[sample(1:nrow(data),size=subsampleSize,replace=TRUE),]
  YC <- replicate_data_set[replicate_data_set$trt==0,]$y
  YT <- replicate_data_set[replicate_data_set$trt==1,]$y
  res <- boundsNoCov_res(0:1, YT, YC, maxBen, maxHarm)
  return(res[1:2])
}

#############################################################################
##Design the true population##
#############################################################################
YCpop <- c(0,1)
YTpop <- c(0,1)
maxBen <- 100
maxHarm <- 0

#############################################################################
##Simulate randomized trials##
#############################################################################
nsim <- 10000
bootrep <- 5000
N <- 1000 ##update
q <- 0.95
k <- 500 ##update
set.seed(9085983)
seed <- sample(10^7,nsim)
nsim <- nsim/10
seed <- seed[4001:5000]

m.candidates <- unique(ceiling(N * q^(0:k)))
m.candidates <- m.candidates[m.candidates>=20]
numCandidates <- length(m.candidates)

bounds <- matrix(data = NA, nrow = nsim, ncol = 3) #bounds and eps
LB.CI_mn <- matrix(data = NA, nrow = nsim, ncol = 3) ##first two elements are limits of CI, third element is the m
UB.CI_mn <- matrix(data = NA, nrow = nsim, ncol = 3)

for (s in 1:nsim){
  set.seed(seed[s])
  YCsamp <- sample(YCpop, size = N/2, replace = TRUE)
  YTsamp <- sample(YTpop, size = N/2, replace = TRUE)
  
  bounds[s,] <- boundsNoCov_res(0:1, YTsamp, YCsamp, maxBen, maxHarm)
  
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

save(bounds, LB.CI_mn, UB.CI_mn, file="N1000_5.Rdata")
