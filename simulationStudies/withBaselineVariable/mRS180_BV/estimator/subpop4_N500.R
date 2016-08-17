##5/31/16

rm(list=ls())

#############################################################################
##Libraries##
#############################################################################
#install.packages("lpSolveAPI")
library(lpSolveAPI)

#############################################################################
##Function for estimating bounds##
#############################################################################
boundsNoCov_res <- function(ordinalScale, YT, YC, maxBen, maxHarm){
  nT <- length(YT) #number of treatment subjects
  nC <- length(YC) #number of control subjects
  
  if(nT == 0 | nC == 0){
    warning("WARNING: YT or YC is empty")
    return(c(0,1,NA))
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
##Marginal distribution of Baseline Variable##
##Marginal distributions of outcome under treatment and control, by stratum##
#############################################################################
setwd("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/mRS180_BV/boundParameters")

load("subpop4.Rdata")

##Convert cdf's to pmf's, since that's easier to use when sampling
pT4.1 <- FT4.1 - c(0,FT4.1[-length(FT4.1)])
pT4.2 <- FT4.2 - c(0,FT4.2[-length(FT4.2)])
pT4.3 <- FT4.3 - c(0,FT4.3[-length(FT4.3)])
pT4.4 <- FT4.4 - c(0,FT4.4[-length(FT4.4)])

pC4.1 <- FC4.1 - c(0,FC4.1[-length(FC4.1)])
pC4.2 <- FC4.2 - c(0,FC4.2[-length(FC4.2)])
pC4.3 <- FC4.3 - c(0,FC4.3[-length(FC4.3)])
pC4.4 <- FC4.4 - c(0,FC4.4[-length(FC4.4)])
#############################################################################
##Simulate randomized trials and compute bound estimates##
#############################################################################
maxBen <- 100
maxHarm <- 100
ordinalScale <- 1:7

nsim <- 10000 
set.seed(2057601)
seed <- sample(10^7,nsim)

N <- 500
K <- 4

bounds <- matrix(data = NA, nrow = nsim, ncol = 3) #bounds and eps, eps should be 0 for all simulated trials

for (m in 1:nsim){
  print(m)
  set.seed(seed[m])
  data <- data.frame(A = c(rep("Medical",N/2), rep("Surgical", N/2)), Y = NA)
  data$X <- sample(1:K, size = nrow(data), replace = TRUE, prob = p4)
  
  data$Y[data$A == "Medical" & data$X == 1] <- sample(ordinalScale, size = sum(data$A == "Medical" & data$X == 1), replace = TRUE, prob = pC4.1)
  data$Y[data$A == "Medical" & data$X == 2] <- sample(ordinalScale, size = sum(data$A == "Medical" & data$X == 2), replace = TRUE, prob = pC4.2)
  data$Y[data$A == "Medical" & data$X == 3] <- sample(ordinalScale, size = sum(data$A == "Medical" & data$X == 3), replace = TRUE, prob = pC4.3)
  data$Y[data$A == "Medical" & data$X == 4] <- sample(ordinalScale, size = sum(data$A == "Medical" & data$X == 4), replace = TRUE, prob = pC4.4)
  
  data$Y[data$A == "Surgical" & data$X == 1] <- sample(ordinalScale, size = sum(data$A == "Surgical" & data$X == 1), replace = TRUE, prob = pT4.1)
  data$Y[data$A == "Surgical" & data$X == 2] <- sample(ordinalScale, size = sum(data$A == "Surgical" & data$X == 2), replace = TRUE, prob = pT4.2)
  data$Y[data$A == "Surgical" & data$X == 3] <- sample(ordinalScale, size = sum(data$A == "Surgical" & data$X == 3), replace = TRUE, prob = pT4.3)
  data$Y[data$A == "Surgical" & data$X == 4] <- sample(ordinalScale, size = sum(data$A == "Surgical" & data$X == 4), replace = TRUE, prob = pT4.4)
  
  bounds1 <- boundsNoCov_res(ordinalScale, data$Y[data$X == 1 & data$A == "Surgical"], data$Y[data$X == 1 & data$A == "Medical"], maxBen, maxHarm) 
  bounds2 <- boundsNoCov_res(ordinalScale, data$Y[data$X == 2 & data$A == "Surgical"], data$Y[data$X == 2 & data$A == "Medical"], maxBen, maxHarm)
  bounds3 <- boundsNoCov_res(ordinalScale, data$Y[data$X == 3 & data$A == "Surgical"], data$Y[data$X == 3 & data$A == "Medical"], maxBen, maxHarm)
  bounds4 <- boundsNoCov_res(ordinalScale, data$Y[data$X == 4 & data$A == "Surgical"], data$Y[data$X == 4 & data$A == "Medical"], maxBen, maxHarm)
  
  bounds[m,] <- mean(data$X == 1) * bounds1 + mean(data$X == 2) * bounds2 + mean(data$X == 3) * bounds3 + mean(data$X == 4) * bounds4
}

bounds <- data.frame(bounds)
names(bounds) <- c("lb","ub","eps")
setwd("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/withBaselineVariable/mRS180_BV/estimator")
save(bounds, file="subpop4_N500.Rdata")
