##5/31/16

rm(list=ls())

#############################################################################
##Libraries##
#############################################################################
#install.packages("lpSolveAPI")
library(lpSolveAPI)

#############################################################################
##Get MISTIE II 180-day mRS data##
#############################################################################
setwd("~/Dropbox/research/data/MISTIEII")
load("MISTIEIIdata.Rdata")

##We only need mRSData and baselineData right now
rm(clotData)

##Also, we want to focus on the data at 180-days
mrs180 <- subset(mrsData, Follow_up_Visit == 180) 
rm(mrsData)
length(unique(mrs180$patientName)) ##96

##Merge mrs180 with baselineData
data <- merge(baselineData, mrs180, by = "patientName")

rm(mrs180, baselineData)

##Only include people with baseline NIHSS and 180-day mRS score
data$patientName[is.na(data$Enrollment_NIHSS_Total)] 
data <- subset(data, !is.na(Enrollment_NIHSS_Total))
nrow(data) ##95

sum(is.na(data$rankin_score)) ##6 people missing rankin_score
data$patientName[is.na(data$rankin_score)] 
data <- subset(data, !is.na(rankin_score)) 
nrow(data) ##89

##We will use the data from these 89 subjects to form the population
##distribution

data <- data.frame(id = data$patientName, 
                   NIHSS = data$Enrollment_NIHSS_Total,
                   rankin_score = data$rankin_score,
                   group = data$Group_Assigned.x)


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

data$rankin_score <- sapply(data$rankin_score, levelfun)

ordinalScale <- 1:7
L <- length(ordinalScale) ##number of levels


#############################################################################
##Marginal distributions for K = 4##
#############################################################################
subpop1 <- subset(data, NIHSS <= 4)
subpop2 <- subset(data, NIHSS > 4 & NIHSS <= 15)
subpop3 <- subset(data, NIHSS > 15 & NIHSS <= 20)
subpop4 <- subset(data, NIHSS > 20)

p4.1 <- nrow(subpop1)/nrow(data)
p4.2 <- nrow(subpop2)/nrow(data)
p4.3 <- nrow(subpop3)/nrow(data)
p4.4 <- nrow(subpop4)/nrow(data)

YT1 <- subpop1$rankin_score[subpop1$group == "Surgical"]
YC1 <- subpop1$rankin_score[subpop1$group == "Medical"]
YT2 <- subpop2$rankin_score[subpop2$group == "Surgical"]
YC2 <- subpop2$rankin_score[subpop2$group == "Medical"]
YT3 <- subpop3$rankin_score[subpop3$group == "Surgical"]
YC3 <- subpop3$rankin_score[subpop3$group == "Medical"]
YT4 <- subpop4$rankin_score[subpop4$group == "Surgical"]
YC4 <- subpop4$rankin_score[subpop4$group == "Medical"]

rm(subpop1, subpop2, subpop3, subpop4)

FC4.1 <- FC4.2 <- FC4.3 <- FC4.4 <- rep(0,L)
for (i in 1:L){
  FC4.1[i] <- sum(YC1 <= ordinalScale[i])/length(YC1)
  FC4.2[i] <- sum(YC2 <= ordinalScale[i])/length(YC2)
  FC4.3[i] <- sum(YC3 <= ordinalScale[i])/length(YC3)
  FC4.4[i] <- sum(YC4 <= ordinalScale[i])/length(YC4)
}

rm(YC1, YC2, YC3, YC4)

FT4.1 <- FT4.2 <- FT4.3 <- FT4.4 <- rep(0,L)
for (i in 1:L){
  FT4.1[i] <- sum(YT1 <= ordinalScale[i])/length(YT1)
  FT4.2[i] <- sum(YT2 <= ordinalScale[i])/length(YT2)
  FT4.3[i] <- sum(YT3 <= ordinalScale[i])/length(YT3)
  FT4.4[i] <- sum(YT4 <= ordinalScale[i])/length(YT4)
}

rm(YT1, YT2, YT3, YT4)

#############################################################################
##Marginal distributions for K = 3##
#############################################################################
p3.1 <- p4.1 + p4.2
p3.2 <- p4.3 
p3.3 <- p4.4


FC3.1 <- p4.1/p3.1*FC4.1 + p4.2/p3.1*FC4.2
FC3.2 <- FC4.3 
FC3.3 <- FC4.4

FT3.1 <- p4.1/p3.1*FT4.1 + p4.2/p3.1*FT4.2
FT3.2 <- FT4.3 
FT3.3 <- FT4.4

#############################################################################
##Marginal distributions for K = 2##
#############################################################################
p2.1 <- p3.1 + p3.2
p2.2 <- p3.3

FC2.1 <- p3.1/p2.1*FC3.1 + p3.2/p2.1*FC3.2
FC2.2 <- FC3.3

FT2.1 <- p3.1/p2.1*FT3.1 + p3.2/p2.1*FT3.2
FT2.2 <- FT3.3



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

bounds1 <- boundsNoCov_res(ordinalScale, FT4.1, FC4.1, maxBen, maxHarm)
bounds2 <- boundsNoCov_res(ordinalScale, FT4.2, FC4.2, maxBen, maxHarm)
bounds3 <- boundsNoCov_res(ordinalScale, FT4.3, FC4.3, maxBen, maxHarm)
bounds4 <- boundsNoCov_res(ordinalScale, FT4.4, FC4.4, maxBen, maxHarm)

res4 <- p4.1*bounds1 + p4.2*bounds2 + p4.3 * bounds3 + p4.4*bounds4

bounds1 <- boundsNoCov_res(ordinalScale, FT3.1, FC3.1, maxBen, maxHarm)
bounds2 <- boundsNoCov_res(ordinalScale, FT3.2, FC3.2, maxBen, maxHarm)
bounds3 <- boundsNoCov_res(ordinalScale, FT3.3, FC3.3, maxBen, maxHarm)

res3 <- p3.1*bounds1 + p3.2*bounds2 + p3.3 * bounds3 

bounds1 <- boundsNoCov_res(ordinalScale, FT2.1, FC2.1, maxBen, maxHarm)
bounds2 <- boundsNoCov_res(ordinalScale, FT2.2, FC2.2, maxBen, maxHarm)

res2 <- p2.1*bounds1 + p2.2*bounds2

#############################################################################
##Save marginal distributions##
#############################################################################
setwd("~/Dropbox/research/github/fraction-who-benefit/simulationStudies/withBaselineVariable/mrs180_BV/boundParameters")

p4 <- c(p4.1, p4.2, p4.3, p4.4)
save(p4, FC4.1, FC4.2, FC4.3, FC4.4, FT4.1, FT4.2, FT4.3, FT4.4, file = "subpop4.Rdata")

p3 <- c(p3.1, p3.2, p3.3)
save(p3, FC3.1, FC3.2, FC3.3, FT3.1, FT3.2, FT3.3, file = "subpop3.Rdata")

p2 <- c(p2.1, p2.2)
save(p2, FC2.1, FC2.2, FT2.1, FT2.2, file = "subpop2.Rdata")

save(res2, res3, res4, file = "trueBounds.Rdata")
