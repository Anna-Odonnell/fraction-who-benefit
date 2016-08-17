##7/12/16

##Compute the distribution of sqrt(n)[psiHAT_l - psi_l(P_n)], under P_n, 
##for n = 100, 1000, 10000, 10^6

########################################################################
##Libraries##
########################################################################
library(lpSolveAPI)

########################################################################
##Function for getting bound estimate##
########################################################################
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
    
    scale <- 1 ##this can be made large to help with solving the linear program
    
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

########################################################################
##Setting inputs##
########################################################################
ordinalScale <- 1:2
maxBen <- 10000
maxHarm <- 10000
nsim <- 10000

########################################################################
##n = 100##
########################################################################
n <- 100
pC1 <- 0.5 + 1/sqrt(n) ##pC1 is P(YC = 1)
pC2 <- 1 - pC1 ##pC2 = P(YC = 2)
pT1 <- 0.5 ##pT1 is P(YT = 1)
pT2 <- 1-pT1 ##pT1 is P(YT = 2)

set.seed(785013)
result100 <- rep(NA, nsim)
for (i in 1:nsim){
  data <- data.frame(A = sample(0:1, size = n, replace = TRUE), Y = NA)
  data$Y[data$A == 0] <- sample(1:2, size = sum(data$A == 0), replace = TRUE, prob = c(pC1,pC2))
  data$Y[data$A == 1] <- sample(1:2, size = sum(data$A == 1), replace = TRUE, prob = c(pT1,pT2))
  result100[i] <- boundsNoCov_res(ordinalScale, data$Y[data$A == 1], data$Y[data$A == 0], maxBen, maxHarm)[1]
}
result100 <- result100*sqrt(n)-1

########################################################################
##n = 1000##
########################################################################
n <- 1000
pC1 <- 0.5 + 1/sqrt(n) ##pC1 is P(YC = 1)
pC2 <- 1 - pC1 ##pC2 = P(YC = 2)
pT1 <- 0.5 ##pT1 is P(YT = 1)
pT2 <- 1-pT1 ##pT1 is P(YT = 2)

set.seed(1165429)
result1000 <- rep(NA, nsim)
for (i in 1:nsim){
  data <- data.frame(A = sample(0:1, size = n, replace = TRUE), Y = NA)
  data$Y[data$A == 0] <- sample(1:2, size = sum(data$A == 0), replace = TRUE, prob = c(pC1,pC2))
  data$Y[data$A == 1] <- sample(1:2, size = sum(data$A == 1), replace = TRUE, prob = c(pT1,pT2))
  result1000[i] <- boundsNoCov_res(ordinalScale, data$Y[data$A == 1], data$Y[data$A == 0], maxBen, maxHarm)[1]
}
result1000 <- result1000*sqrt(n)-1

########################################################################
##n = 10000##
########################################################################
n <- 10000
pC1 <- 0.5 + 1/sqrt(n) ##pC1 is P(YC = 1)
pC2 <- 1 - pC1 ##pC2 = P(YC = 2)
pT1 <- 0.5 ##pT1 is P(YT = 1)
pT2 <- 1-pT1 ##pT1 is P(YT = 2)

set.seed(7692624)
result10000 <- rep(NA, nsim)
for (i in 1:nsim){
  data <- data.frame(A = sample(0:1, size = n, replace = TRUE), Y = NA)
  data$Y[data$A == 0] <- sample(1:2, size = sum(data$A == 0), replace = TRUE, prob = c(pC1,pC2))
  data$Y[data$A == 1] <- sample(1:2, size = sum(data$A == 1), replace = TRUE, prob = c(pT1,pT2))
  result10000[i] <- boundsNoCov_res(ordinalScale, data$Y[data$A == 1], data$Y[data$A == 0], maxBen, maxHarm)[1]
}
result10000 <- result10000*sqrt(n)-1

########################################################################
##n = 10^6##
########################################################################
n <- 10^6
pC1 <- 0.5 + 1/sqrt(n) ##pC1 is P(YC = 1)
pC2 <- 1 - pC1 ##pC2 = P(YC = 2)
pT1 <- 0.5 ##pT1 is P(YT = 1)
pT2 <- 1-pT1 ##pT1 is P(YT = 2)

set.seed(9816388)
resultMillion <- rep(NA, nsim)
for (i in 1:nsim){
  data <- data.frame(A = sample(0:1, size = n, replace = TRUE), Y = NA)
  data$Y[data$A == 0] <- sample(1:2, size = sum(data$A == 0), replace = TRUE, prob = c(pC1,pC2))
  data$Y[data$A == 1] <- sample(1:2, size = sum(data$A == 1), replace = TRUE, prob = c(pT1,pT2))
  resultMillion[i] <- boundsNoCov_res(ordinalScale, data$Y[data$A == 1], data$Y[data$A == 0], maxBen, maxHarm)[1]
}
resultMillion <- resultMillion*sqrt(n)-1

########################################################################
##Save results##
########################################################################
save(result100, result1000, result10000, resultMillion, file = "distributions_Pn.Rdata")
