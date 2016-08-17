##R code

library(lpSolveAPI)


boundsNoCov_res <- function(ordinalScale, YT, YC, maxBen, maxHarm){
  nT <- length(YT) #number of treatment subjects
  nC <- length(YC) #number of control subjects
  
  if(nT == 0 | nC == 0){
    stop("ERROR: YT or YC is empty")
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



boundsCov_res <- function(ordinalScale, YT1, YC1, YT2, YC2, maxBen, maxHarm){
  x1 <- boundsNoCov_res(ordinalScale, YT1, YC1, maxBen, maxHarm)
  x2 <- boundsNoCov_res(ordinalScale, YT2, YC2, maxBen, maxHarm)
  p1 <- (length(YT1) + length(YC1))/(length(YT1)+length(YC1)+length(YT2)+length(YC2))
  lb <- x1[1]*p1 + x2[1]*(1-p1)
  ub <- x1[2]*p1 + x2[2]*(1-p1)
  return(c(lb,ub,x1,x2))
}

##The code above was used for our analyses because support restrictions were
##applied to the entire population. If your support restrictions differ between the two subpopulations,
##please use the code below.

# boundsCov_res <- function(ordinalScale, YT1, YC1, YT2, YC2, maxBen1, maxHarm1, maxBen2, maxHarm2){
#   x1 <- boundsNoCov_res(ordinalScale, YT1, YC1, maxBen1, maxHarm1)
#   x2 <- boundsNoCov_res(ordinalScale, YT2, YC2, maxBen2, maxHarm2)
#   p1 <- (length(YT1) + length(YC1))/(length(YT1)+length(YC1)+length(YT2)+length(YC2))
#   lb <- x1[1]*p1 + x2[1]*(1-p1)
#   ub <- x1[2]*p1 + x2[2]*(1-p1)
#   return(c(lb,ub,x1,x2))
# }
