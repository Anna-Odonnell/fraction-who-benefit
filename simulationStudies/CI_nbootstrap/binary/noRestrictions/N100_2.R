rm(list=ls())

#############################################################################
##Libraries##
#############################################################################
library(boot)
library(lpSolveAPI)

#############################################################################
##Function for estimating bounds##
#############################################################################
boundsNoCov_res <- function(ordinalScale, YT, YC, maxBen, maxHarm){
  nT <- length(YT) #number of treatment subjects
  nC <- length(YC) #number of control subjects
  
  if(nT == 0 | nC == 0){
    stop("WARNING: YT or YC is empty")
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
##Function for n-bootstrap##
#############################################################################

compute_estimators_on_bootstrap_replicate_of_dataset <- function(data,index,maxBen,maxHarm) {
  replicate_data_set <- data[index,]
  YC <- replicate_data_set[replicate_data_set$trt==0,]$y
  YT <- replicate_data_set[replicate_data_set$trt==1,]$y
  res <- boundsNoCov_res(0:1, YT, YC, maxBen, maxHarm)
  return(res)
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
nsim <- 10000 
set.seed(4983095)
seed <- sample(10^7, nsim)
seed <- seed[2001:4000]
nsim <- nsim/5
bootrep <- 10000
N <- 100 

bounds <- matrix(data = NA, nrow = nsim, ncol = 3) #bounds and eps
LB.CI_perc <- matrix(data = NA, nrow = nsim, ncol = 2)
UB.CI_perc <- matrix(data = NA, nrow = nsim, ncol = 2)

for (m in 1:nsim){
  set.seed(seed[m])
  YCsamp <- sample(YCpop, size = N/2, replace = TRUE)
  YTsamp <- sample(YTpop, size = N/2, replace = TRUE)
  datasamp <- data.frame(trt=c(rep(1,length(YTsamp)),rep(0,length(YCsamp))),y=c(YTsamp,YCsamp)) #trt is 1 if treatment and 0 if control, y is the outcome
  
  boot.out <- boot(data = datasamp, statistic= compute_estimators_on_bootstrap_replicate_of_dataset, R=bootrep, maxBen = maxBen, maxHarm = maxHarm)
  bounds[m,] <- boot.out$t0

  LB.CI_perc[m,] <- boot.ci(boot.out, index = 1, conf=0.95, type="perc")$perc[4:5] 
  UB.CI_perc[m,] <- boot.ci(boot.out, index = 2, conf=0.95, type="perc")$perc[4:5]
}

save(bounds, LB.CI_perc, UB.CI_perc, file="N100_2.Rdata")
