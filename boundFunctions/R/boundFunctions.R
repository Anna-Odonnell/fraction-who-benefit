##R code

library(lpSolveAPI)


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



boundsCov_res <- function(YT1, YC1, YT2, YC2, maxBen, maxHarm){
  x1 <- boundsNoCov_res(YT1, YC1, maxBen, maxHarm)
  x2 <- boundsNoCov_res(YT2, YC2, maxBen, maxHarm)
  p1 <- (length(YT1) + length(YC1))/(length(YT1)+length(YC1)+length(YT2)+length(YC2))
  lb <- x1[1]*p1 + x2[1]*(1-p1)
  ub <- x1[2]*p1 + x2[2]*(1-p1)
  return(c(lb,ub,x1,x2))
}

##The code above was convenient for our analyses because all support restrictions
##applied to the entire population. If your support restrictions differ between the two subpopulations,
##please use the code below.

# boundsCov_res <- function(YT1, YC1, YT2, YC2, maxBen1, maxHarm1, maxBen2, maxHarm2){
#   x1 <- boundsNoCov_res(YT1, YC1, maxBen1, maxHarm1)
#   x2 <- boundsNoCov_res(YT2, YC2, maxBen2, maxHarm2)
#   p1 <- (length(YT1) + length(YC1))/(length(YT1)+length(YC1)+length(YT2)+length(YC2))
#   lb <- x1[1]*p1 + x2[1]*(1-p1)
#   ub <- x1[2]*p1 + x2[2]*(1-p1)
#   return(c(lb,ub,x1,x2))
# }
