##5/23/16

rm(list=ls())

#############################################################################
##Libraries##
#############################################################################
#install.packages("lpSolveAPI")
library(lpSolveAPI)

#############################################################################
##Gather the RICV5 data##
#############################################################################
setwd("~/Dropbox/research/data/MISTIEII")
load("MISTIEIIdata.Rdata")

##We can ignore mrsData
rm(mrsData)

##Merge baselineData and clotData
data <- merge(baselineData, clotData, by = "patientName")

rm(baselineData, clotData)

##Calculate reduction in clot volume
data$volchange <- data$Pre_Rand_ICH_Volume_RC - data$eot_ich_9_13

##Focus on the variables of interest
data <- data.frame(data$patientName, data$Group_Assigned, data$volchange, data$Pre_Rand_ICH_Volume_RC)
names(data) <- c("id", "group", "RICV", "baseVol")

##Discretize RICV using a bin length of 5 mL##
levelfun <- function(x){ 
  if(x < 0){
    y = 1
  } else if(x >= 0 & x < 5){
    y = 2
  } else if(x >= 5 & x < 10){
    y = 3
  } else if(x >= 10 & x < 15){
    y = 4
  } else if(x >= 15 & x < 20){
    y = 5
  } else {
    y = 6
  } 
  return(y)
}

data$RICV5 <- sapply(data$RICV, levelfun)

rm(levelfun)

ordinalScale <- 1:6
L <- length(ordinalScale)

#############################################################################
##Marginal distributions for K = 8##
#############################################################################
q1 <- quantile(data$baseVol, probs = 1/8)
q2 <- quantile(data$baseVol, probs = 2/8)
q3 <- quantile(data$baseVol, probs = 3/8)
q4 <- quantile(data$baseVol, probs = 4/8)
q5 <- quantile(data$baseVol, probs = 5/8)
q6 <- quantile(data$baseVol, probs = 6/8)
q7 <- quantile(data$baseVol, probs = 7/8)

subpop1 <- subset(data, baseVol <= q1)
subpop2 <- subset(data, baseVol > q1 & baseVol <= q2)
subpop3 <- subset(data, baseVol > q2 & baseVol <= q3)
subpop4 <- subset(data, baseVol > q3 & baseVol <= q4)
subpop5 <- subset(data, baseVol > q4 & baseVol <= q5)
subpop6 <- subset(data, baseVol > q5 & baseVol <= q6)
subpop7 <- subset(data, baseVol > q6 & baseVol <= q7)
subpop8 <- subset(data, baseVol > q7)

rm(q1, q2, q3, q4, q5, q6, q7)

p8.1 <- nrow(subpop1)/nrow(data)
p8.2 <- nrow(subpop2)/nrow(data)
p8.3 <- nrow(subpop3)/nrow(data)
p8.4 <- nrow(subpop4)/nrow(data)
p8.5 <- nrow(subpop5)/nrow(data)
p8.6 <- nrow(subpop6)/nrow(data)
p8.7 <- nrow(subpop7)/nrow(data)
p8.8 <- nrow(subpop8)/nrow(data)

YT1 <- subpop1$RICV5[subpop1$group == "Surgical"]
YC1 <- subpop1$RICV5[subpop1$group == "Medical"]
YT2 <- subpop2$RICV5[subpop2$group == "Surgical"]
YC2 <- subpop2$RICV5[subpop2$group == "Medical"]
YT3 <- subpop3$RICV5[subpop3$group == "Surgical"]
YC3 <- subpop3$RICV5[subpop3$group == "Medical"]
YT4 <- subpop4$RICV5[subpop4$group == "Surgical"]
YC4 <- subpop4$RICV5[subpop4$group == "Medical"]
YT5 <- subpop5$RICV5[subpop5$group == "Surgical"]
YC5 <- subpop5$RICV5[subpop5$group == "Medical"]
YT6 <- subpop6$RICV5[subpop6$group == "Surgical"]
YC6 <- subpop6$RICV5[subpop6$group == "Medical"]
YT7 <- subpop7$RICV5[subpop7$group == "Surgical"]
YC7 <- subpop7$RICV5[subpop7$group == "Medical"]
YT8 <- subpop8$RICV5[subpop8$group == "Surgical"]
YC8 <- subpop8$RICV5[subpop8$group == "Medical"]

rm(subpop1, subpop2, subpop3, subpop4, subpop5, subpop6, subpop7, subpop8)

FC8.1 <- FC8.2 <- FC8.3 <- FC8.4 <- FC8.5 <- FC8.6 <- FC8.7 <- FC8.8 <- rep(0,L)
for (i in 1:L){
  FC8.1[i] <- sum(YC1 <= ordinalScale[i])/length(YC1)
  FC8.2[i] <- sum(YC2 <= ordinalScale[i])/length(YC2)
  FC8.3[i] <- sum(YC3 <= ordinalScale[i])/length(YC3)
  FC8.4[i] <- sum(YC4 <= ordinalScale[i])/length(YC4)
  FC8.5[i] <- sum(YC5 <= ordinalScale[i])/length(YC5)
  FC8.6[i] <- sum(YC6 <= ordinalScale[i])/length(YC6)
  FC8.7[i] <- sum(YC7 <= ordinalScale[i])/length(YC7)
  FC8.8[i] <- sum(YC8 <= ordinalScale[i])/length(YC8)
}

rm(YC1, YC2, YC3, YC4, YC5, YC6, YC7, YC8)

FT8.1 <- FT8.2 <- FT8.3 <- FT8.4 <- FT8.5 <- FT8.6 <- FT8.7 <- FT8.8 <- rep(0,L)
for (i in 1:L){
  FT8.1[i] <- sum(YT1 <= ordinalScale[i])/length(YT1)
  FT8.2[i] <- sum(YT2 <= ordinalScale[i])/length(YT2)
  FT8.3[i] <- sum(YT3 <= ordinalScale[i])/length(YT3)
  FT8.4[i] <- sum(YT4 <= ordinalScale[i])/length(YT4)
  FT8.5[i] <- sum(YT5 <= ordinalScale[i])/length(YT5)
  FT8.6[i] <- sum(YT6 <= ordinalScale[i])/length(YT6)
  FT8.7[i] <- sum(YT7 <= ordinalScale[i])/length(YT7)
  FT8.8[i] <- sum(YT8 <= ordinalScale[i])/length(YT8)
}

rm(YT1, YT2, YT3, YT4, YT5, YT6, YT7, YT8)

#############################################################################
##Marginal distributions for K = 4##
#############################################################################
p4.1 <- p8.1 + p8.2
p4.2 <- p8.3 + p8.4
p4.3 <- p8.5 + p8.6
p4.4 <- p8.7 + p8.8

FC4.1 <- p8.1/p4.1*FC8.1 + p8.2/p4.1*FC8.2
FC4.2 <- p8.3/p4.2*FC8.3 + p8.4/p4.2*FC8.4
FC4.3 <- p8.5/p4.3*FC8.5 + p8.6/p4.3*FC8.6
FC4.4 <- p8.7/p4.4*FC8.7 + p8.8/p4.4*FC8.8

FT4.1 <- p8.1/p4.1*FT8.1 + p8.2/p4.1*FT8.2
FT4.2 <- p8.3/p4.2*FT8.3 + p8.4/p4.2*FT8.4
FT4.3 <- p8.5/p4.3*FT8.5 + p8.6/p4.3*FT8.6
FT4.4 <- p8.7/p4.4*FT8.7 + p8.8/p4.4*FT8.8

#############################################################################
##Marginal distributions for K = 2##
#############################################################################
p2.1 <- p4.1 + p4.2
p2.2 <- p4.3 + p4.4

FC2.1 <- p4.1/p2.1*FC4.1 + p4.2/p2.1*FC4.2
FC2.2 <- p4.3/p2.2*FC4.3 + p4.4/p2.2*FC4.4

FT2.1 <- p4.1/p2.1*FT4.1 + p4.2/p2.1*FT4.2
FT2.2 <- p4.3/p2.2*FT4.3 + p4.4/p2.2*FT4.4




#############################################################################
##Function for computing bounds##
#############################################################################
boundsNoCov_res <- function(ordinalScale, cdf_T, cdf_C, maxBen, maxHarm){
  ordinalScale <- sort(ordinalScale, decreasing = FALSE)
  L <- length(ordinalScale)
  varCount <- L^2 #number of pi_{i,j}'s
  
  #The matrix of pi_{i,j}'s is a L x L matrix with varCount pi_{i,j}'s.
  scale <- 1
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


#############################################################################
##Getting bound parameters##
#############################################################################
maxBen <- 100
maxHarm <- 100

bounds1 <- boundsNoCov_res(ordinalScale, FT8.1, FC8.1, maxBen, maxHarm)
bounds2 <- boundsNoCov_res(ordinalScale, FT8.2, FC8.2, maxBen, maxHarm)
bounds3 <- boundsNoCov_res(ordinalScale, FT8.3, FC8.3, maxBen, maxHarm)
bounds4 <- boundsNoCov_res(ordinalScale, FT8.4, FC8.4, maxBen, maxHarm)
bounds5 <- boundsNoCov_res(ordinalScale, FT8.5, FC8.5, maxBen, maxHarm)
bounds6 <- boundsNoCov_res(ordinalScale, FT8.6, FC8.6, maxBen, maxHarm)
bounds7 <- boundsNoCov_res(ordinalScale, FT8.7, FC8.7, maxBen, maxHarm)
bounds8 <- boundsNoCov_res(ordinalScale, FT8.8, FC8.8, maxBen, maxHarm)

res8 <- p8.1*bounds1 + p8.2*bounds2 + p8.3 * bounds3 + p8.4*bounds4 + p8.5*bounds5 + p8.6*bounds6 + p8.7*bounds7 + p8.8*bounds8

bounds1 <- boundsNoCov_res(ordinalScale, FT4.1, FC4.1, maxBen, maxHarm)
bounds2 <- boundsNoCov_res(ordinalScale, FT4.2, FC4.2, maxBen, maxHarm)
bounds3 <- boundsNoCov_res(ordinalScale, FT4.3, FC4.3, maxBen, maxHarm)
bounds4 <- boundsNoCov_res(ordinalScale, FT4.4, FC4.4, maxBen, maxHarm)

res4 <- p4.1*bounds1 + p4.2*bounds2 + p4.3 * bounds3 + p4.4*bounds4

bounds1 <- boundsNoCov_res(ordinalScale, FT2.1, FC2.1, maxBen, maxHarm)
bounds2 <- boundsNoCov_res(ordinalScale, FT2.2, FC2.2, maxBen, maxHarm)

res2 <- p2.1*bounds1 + p2.2*bounds2

#############################################################################
##Save marginal distributions##
#############################################################################
setwd("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/withBaselineVariable/RICV5_BV/boundParameters")

p8 <- c(p8.1, p8.2, p8.3, p8.4, p8.5, p8.6, p8.7, p8.8)
save(p8, FC8.1, FC8.2, FC8.3, FC8.4, FC8.5, FC8.6, FC8.7, FC8.8, FT8.1, FT8.2, FT8.3, FT8.4, FT8.5, FT8.6, FT8.7, FT8.8, file = "subpop8.Rdata")

p4 <- c(p4.1, p4.2, p4.3, p4.4)
save(p4, FC4.1, FC4.2, FC4.3, FC4.4, FT4.1, FT4.2, FT4.3, FT4.4, file = "subpop4.Rdata")

p2 <- c(p2.1, p2.2)
save(p2, FC2.1, FC2.2, FT2.1, FT2.2, file = "subpop2.Rdata")

save(res2, res4, res8, file = "trueBounds.Rdata")
