##This code is for making Tables 1a and 1b of the paper. We
##analyze the distribution of the lower and upper bound estimators and
##the performance of the CI methods (n and m-out-of-n). 

rm(list=ls())

library(lpSolveAPI)
library(xtable)
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
##Empty tables##
#############################################################################

tab1 <- matrix(NA, nrow = 6, ncol = 6)
tab2 <- matrix(NA, nrow = 6, ncol = 12)


#############################################################################
##Reduction in clot volume##
#############################################################################
##First calculate the lower and upper bound parameters
##YCd.txt and YTd.txt were created in changeInVol.R in the MISTIEanalysis folder
YCpop <- unlist(read.table("YCd.txt"))
YTpop <- unlist(read.table("YTd.txt"))
maxBen <- 100
maxHarm <- 100
temp <- boundsNoCov_res(YTpop, YCpop, maxBen, maxHarm)
trueLB <- temp[1]
trueUB <- temp[2]
rm(temp, maxBen, maxHarm, YCpop, YTpop, boundsNoCov_res)

##Go to CI_nbootstrap/RICV folder

##LOAD "results_clotRed_noRes_N100_Perc.Rdata"
LB <- bounds[,1]
UB <- bounds[,2]
rm(bounds)
tab1[1,1] <- mean(LB)-trueLB ##Bias
tab1[1,2] <- mean(UB)-trueUB
tab1[4,1] <- sd(LB) ##SE
tab1[4,2] <- sd(UB)
tab2[1,1] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[1,3] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[4,1] <- mean(LB.CI[,2] - LB.CI[,1])
tab2[4,3] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB,UB,LB.CI,UB.CI)

##LOAD "results_clotRed_noRes_N500_Perc.Rdata"
LB <- bounds[,1]
UB <- bounds[,2]
rm(bounds)
tab1[2,1] <- mean(LB)-trueLB ##Bias
tab1[2,2] <- mean(UB)-trueUB
tab1[5,1] <- sd(LB) ##SE
tab1[5,2] <- sd(UB)
tab2[2,1] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[2,3] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[5,1] <- mean(LB.CI[,2] - LB.CI[,1]) #Average width
tab2[5,3] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB,UB,LB.CI,UB.CI)

##LOAD "results_clotRed_noRes_N1000_Perc.Rdata"
LB <- bounds[,1]
UB <- bounds[,2]
rm(bounds)
tab1[3,1] <- mean(LB)-trueLB ##Bias
tab1[3,2] <- mean(UB)-trueUB
tab1[6,1] <- sd(LB) ##SE
tab1[6,2] <- sd(UB)
tab2[3,1] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[3,3] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[6,1] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[6,3] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB,UB,LB.CI,UB.CI)

##Go to CI_mofn/RICV folder

##LOAD "results_clotVol_noRes_N100_mofn.Rdata"
LB.CI <- LB.CI_mn
UB.CI <- UB.CI_mn
tab2[1,2] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[1,4] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[4,2] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[4,4] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB.CI,UB.CI,LB.CI_mn,UB.CI_mn,bounds)

##LOAD "results_clotVol_noRes_N500_mofn.Rdata"
LB.CI <- LB.CI_mn
UB.CI <- UB.CI_mn
tab2[2,2] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[2,4] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[5,2] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[5,4] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB.CI,UB.CI,LB.CI_mn,UB.CI_mn,bounds)

##LOAD "results_clotVol_noRes_N1000_mofn.Rdata"
LB.CI <- LB.CI_mn
UB.CI <- UB.CI_mn
tab2[3,2] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[3,4] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[6,2] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[6,4] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB.CI,UB.CI,LB.CI_mn,UB.CI_mn,bounds)

#############################################################################
##Binary (no restrictions)##
#############################################################################
trueLB <- 0
trueUB <- 0.5

##Go to CI_nbootstrap/binary/noRestrictions folder

##LOAD "results_binary_noRes_N100_Perc.Rdata"
LB <- bounds[,1]
UB <- bounds[,2]
rm(bounds)
tab1[1,3] <- mean(LB)-trueLB ##Bias
tab1[1,4] <- mean(UB)-trueUB
tab1[4,3] <- sd(LB) ##SE
tab1[4,4] <- sd(UB)
tab2[1,5] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[1,7] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[4,5] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[4,7] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB,UB,LB.CI,UB.CI)

##LOAD "results_binary_noRes_N500_Perc.Rdata"
LB <- bounds[,1]
UB <- bounds[,2]
rm(bounds)
tab1[2,3] <- mean(LB)-trueLB ##Bias
tab1[2,4] <- mean(UB)-trueUB
tab1[5,3] <- sd(LB) ##SE
tab1[5,4] <- sd(UB)
tab2[2,5] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[2,7] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[5,5] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[5,7] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB,UB,LB.CI,UB.CI)

##LOAD "results_binary_noRes_N1000_Perc.Rdata"
LB <- bounds[,1]
UB <- bounds[,2]
rm(bounds)
tab1[3,3] <- mean(LB)-trueLB ##Bias
tab1[3,4] <- mean(UB)-trueUB
tab1[6,3] <- sd(LB) ##SE
tab1[6,4] <- sd(UB)
tab2[3,5] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[3,7] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[6,5] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[6,7] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB,UB,LB.CI,UB.CI)

##Go to CI_mofn/binary/noRestrictions folder

##LOAD "results_binary_noRes_N100_mofn.Rdata"
tab2[1,6] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[1,8] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[4,6] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[4,8] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB.CI,UB.CI,bounds)

##LOAD "results_binary_noRes_N500_mofn.Rdata"
tab2[2,6] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[2,8] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[5,6] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[5,8] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB.CI,UB.CI,bounds)

##LOAD "results_binary_noRes_N1000_mofn.Rdata"
tab2[3,6] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[3,8] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[6,6] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[6,8] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB.CI,UB.CI,bounds)


#############################################################################
##Binary (no harm)##
#############################################################################
trueLB <- 0
trueUB <- 0

##Go to CI_nbootstrap/binary/noHarm folder

##LOAD "results_binary_noHarm_N100_Perc.Rdata"
LB <- bounds[,1]
UB <- bounds[,2]
rm(bounds)
tab1[1,5] <- mean(LB)-trueLB ##Bias
tab1[1,6] <- mean(UB)-trueUB
tab1[4,5] <- sd(LB) ##SE
tab1[4,6] <- sd(UB)
tab2[1,9] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[1,11] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[4,9] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[4,11] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB,UB,LB.CI,UB.CI)

##LOAD "results_binary_noHarm_N500_Perc.Rdata"
LB <- bounds[,1]
UB <- bounds[,2]
rm(bounds)
tab1[2,5] <- mean(LB)-trueLB ##Bias
tab1[2,6] <- mean(UB)-trueUB
tab1[5,5] <- sd(LB) ##SE
tab1[5,6] <- sd(UB)
tab2[2,9] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[2,11] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[5,9] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[5,11] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB,UB,LB.CI,UB.CI)

##LOAD "results_binary_noHarm_N1000_Perc.Rdata"
LB <- bounds[,1]
UB <- bounds[,2]
rm(bounds)
tab1[3,5] <- mean(LB)-trueLB ##Bias
tab1[3,6] <- mean(UB)-trueUB
tab1[6,5] <- sd(LB) ##SE
tab1[6,6] <- sd(UB)
tab2[3,9] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[3,11] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[6,9] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[6,11] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB,UB,LB.CI,UB.CI)

##Go to CI_mofn/binary/noHarm folder

##LOAD "results_binary_noHarm_N100_mofn.Rdata"
LB.CI <- LB.CI_mn
UB.CI <- UB.CI_mn
tab2[1,10] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[1,12] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[4,10] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[4,12] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB.CI,UB.CI,LB.CI_mn,UB.CI_mn,bounds)

##LOAD "results_binary_noHarm_N500_mofn.Rdata"
LB.CI <- LB.CI_mn
UB.CI <- UB.CI_mn
tab2[2,10] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[2,12] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[5,10] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[5,12] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB.CI,UB.CI,LB.CI_mn,UB.CI_mn,bounds)

##LOAD "results_binary_noHarm_N1000_mofn.Rdata")
LB.CI <- LB.CI_mn
UB.CI <- UB.CI_mn
tab2[3,10] <- mean(trueLB >= LB.CI[,1] & trueLB <= LB.CI[,2]) ##Coverage
tab2[3,12] <- mean(trueUB >= UB.CI[,1] & trueUB <= UB.CI[,2])
tab2[6,10] <- mean(LB.CI[,2] - LB.CI[,1]) ##Average width
tab2[6,12] <- mean(UB.CI[,2] - UB.CI[,1])
rm(LB.CI,UB.CI,LB.CI_mn,UB.CI_mn,bounds)

##Completed tables
tab1
tab2

##Get the tables in Latex format
tab1 <- data.frame(tab1)
names(tab1) <- c("RCV-LB", "RCV-UB", "BinaryNoRes-LB", "BinaryNoRes-UB", "BinaryNoHarm-LB", "BinaryNoHarm-UB")
row.names(tab1) <- c("Bias-n=100","Bias-n=500","Bias-n=1000","SE-n=100","SE-n=500","SE-n=1000")
table1 <- xtable(tab1, align = "ccccccc")
digits(table1) <- 3
print(table1)

tab2 <- data.frame(tab2)
names(tab2) <- c("RCV-LB-n","RCV-LB-m","RCV-UB-n","RCV-UB-m","BinaryNoRes-LB-n","BinaryNoRes-LB-m","BinaryNoRes-UB-n","BinaryNoRes-UB-m","BinaryNoHarm-LB-n","BinaryNoHarm-LB-m","BinaryNoHarm-UB-n","BinaryNoHarm-UB-m")
row.names(tab2) <- c("Coverage-n=100","Coverage-n=500","Coverage-n=1000","AvgWid-n=100","AvgWid-n=500","AvgWid-n=1000")
table2 <- xtable(tab2, align = "ccccccccccccc")
digits(table2) <- 3
print(table2)
