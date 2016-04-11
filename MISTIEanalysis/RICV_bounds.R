##Analysis below includes treatment and control patients from MISTIE II, and control patients from ICES group.
##Patients are included only if they had random treatment assignment, non-missing reduction in clot volume, 
##and non-missing baseline clot volume

rm(list=ls())

##########################################################################################################
##Gather the reduction in clot volume data##
##########################################################################################################

load("MISTIEIIdata.Rdata")
##Each of these datasets (i.e., mRSData, clotData, baselineData) has 96 people (42 medical, 54 surgical)
##The people in these datasets meet all inclusion criteria, except having non-missing reduction in clot volume 
##and baseline clot volume

##We can ignore mrsData
rm(mrsData)

##Merge baselineData and clotData
data <- merge(baselineData, clotData, by = "patientName")

rm(baselineData, clotData)

##Check that merge was successful
all(data$Group_Assigned == data$group_assigned) ##true

##Everyone has baseline and EOT reduction in clot volume
sum(is.na(data$Pre_Rand_ICH_Volume_RC)) ##0
sum(is.na(data$eot_ich_9_13)) ##0

##Calculate reduction in clot volume
data$volchange <- data$Pre_Rand_ICH_Volume_RC - data$eot_ich_9_13

##There will be 96 patients in the RICV analysis.

##Focus on the variables of interest
data <- data.frame(data$patientName, data$Group_Assigned, data$volchange, data$Pre_Rand_ICH_Volume_RC)
names(data) <- c("id", "group", "volChange", "baseVol")

control <- subset(data, group == "Medical") ##n = 42
tmt <- subset(data, group == "Surgical") ##n = 54

YT <- tmt$volChange
YC <- control$volChange


m <- median(data$baseVol) ##median baseline clot volume is 43.204999925

control1 <- subset(control, baseVol < m) ##n = 22
tmt1 <- subset(tmt, baseVol < m) ##n = 26
control2 <- subset(control, baseVol >= m) ##n = 20
tmt2 <- subset(tmt, baseVol >= m) ##n = 28

YT1 <- tmt1$volChange
YT2 <- tmt2$volChange
YC1 <- control1$volChange
YC2 <- control2$volChange

##########################################################################################################
##Discretize RICV using a bin length of 5 mL and get bound estimates##
##########################################################################################################
##convert continuous clot volume into a discrete variable
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

YCd <- sapply(YC, levelfun)
YTd <- sapply(YT, levelfun)
YC1d <- sapply(YC1, levelfun)
YT1d <- sapply(YT1, levelfun)
YC2d <- sapply(YC2, levelfun)
YT2d <- sapply(YT2, levelfun)

##Get the bound estimates
ordinalScale <- 1:6

result.noBV <- matrix(0,nrow=10,ncol=3)
result.BV <- matrix(NA,nrow=10,ncol=8)
maxBen <- c(100, 4, 3, 2, 1, 100, 100, 100, 100, 100)
maxHarm <- c(100, 100, 100, 100, 100, 4, 3, 2, 1, 0)

for (i in 1:10){
  result.noBV[i,] <- boundsNoCov_res(ordinalScale, YTd, YCd, maxBen[i], maxHarm[i])
  result.BV[i,] <- boundsCov_res(ordinalScale, YT1d, YC1d, YT2d, YC2d, maxBen[i], maxHarm[i])
}

result.pop <- cbind(result.noBV,result.BV[,1:2])
result.subpop <- result.BV[,3:8]

##The tables below are presented in the Supplementary Materials
result.pop <- data.frame(result.pop)
names(result.pop) <- c("l", "u", "epsilon", "l (using BV)", "u (using BV)")
row.names(result.pop) <- c("No assumptions", "Benefit at most 4", "at most 3", "at most 2", "at most 1", 
                        "Harm at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(result.pop, align = "rccccc")
digits(tab) <- 2
print(tab)

result.subpop <- data.frame(result.subpop)
names(result.subpop) <- c("l", "u", "epsilon", "l", "u", "epsilon")
row.names(result.subpop) <- c("No assumptions", "Benefit at most 4", "at most 3", "at most 2", "at most 1", 
                           "Harm at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(result.subpop, align = "rcccccc")
digits(tab) <- 2
print(tab)

##########################################################################################################
##Discretize RICV using a bin length of 2 mL and get bound estimates##
##########################################################################################################
levelfun <- function(x){ 
  if(x < 0){
    y = 1
  } else if(x >= 0 & x < 2){
    y = 2
  } else if(x >= 2 & x < 4){
    y = 3
  } else if(x >= 4 & x < 6){
    y = 4
  } else if(x >= 6 & x < 8){
    y = 5
  } else if(x >= 8 & x < 10){
    y = 6
  } else if(x >= 10 & x < 12){
    y = 7
  } else if(x >= 12 & x < 14){
    y = 8
  } else if(x >= 14 & x < 16){
    y = 9
  } else if(x >= 16 & x < 18){
    y = 10
  } else if(x >= 18 & x < 20){
    y = 11
  } else {
    y = 12
  } 
  return(y)
}

YCd <- sapply(YC, levelfun)
YTd <- sapply(YT, levelfun)
YC1d <- sapply(YC1, levelfun)
YT1d <- sapply(YT1, levelfun)
YC2d <- sapply(YC2, levelfun)
YT2d <- sapply(YT2, levelfun)

ordinalScale <- 1:12

result.noBV <- matrix(NA,nrow=22,ncol=3)
result.BV <- matrix(NA,nrow=22,ncol=8)
maxBen <- c(100, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
maxHarm <- c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

for (i in 1:22){
  result.noBV[i,] <- boundsNoCov_res(ordinalScale, YTd, YCd, maxBen[i], maxHarm[i])
  result.BV[i,] <- boundsCov_res(ordinalScale, YT1d, YC1d, YT2d, YC2d, maxBen[i], maxHarm[i])
}

result.pop <- cbind(result.noBV,result.BV[,1:2])
result.subpop <- result.BV[,3:8]

##The tables below are presented in the Supplementary Materials
result.pop <- data.frame(result.pop)
names(result.pop) <- c("l", "u", "epsilon", "l (using BV)", "u (using BV)")
row.names(result.pop) <- c("No assumptions", "Benefit at most 10", "at most 9", "at most 8", "at most 7",
                           "at most 6", "at most 5", "at most 4", "at most 3", "at most 2", "at most 1",
                           "Harm at most 10", " at most 9", " at most 8", " at most 7", " at most 6", 
                           " at most 5", " at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(result.pop, align = "rccccc")
digits(tab) <- 2
print(tab)

result.subpop <- data.frame(result.subpop)
names(result.subpop) <- c("l", "u", "epsilon", "l", "u", "epsilon")
row.names(result.subpop) <- c("No assumptions", "Benefit at most 10", "at most 9", "at most 8", "at most 7",
                              "at most 6", "at most 5", "at most 4", "at most 3", "at most 2", "at most 1",
                              "Harm at most 10", " at most 9", " at most 8", " at most 7", " at most 6", 
                              " at most 5", " at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(result.subpop, align = "rcccccc")
digits(tab) <- 2
print(tab)
##########################################################################################################
##Discretize RICV using a bin length of 10 mL and get bound estimates##
##########################################################################################################
levelfun <- function(x){ 
  if(x < 0){
    y = 1
  } else if(x >= 0 & x < 10){
    y = 2
  } else if(x >= 10 & x < 20){
    y = 3
  } else {
    y = 4
  } 
  return(y)
}

YCd <- sapply(YC, levelfun)
YTd <- sapply(YT, levelfun)
YC1d <- sapply(YC1, levelfun)
YT1d <- sapply(YT1, levelfun)
YC2d <- sapply(YC2, levelfun)
YT2d <- sapply(YT2, levelfun)

ordinalScale <- 1:4

result.noBV <- matrix(0,nrow=6,ncol=3)
result.BV <- matrix(NA,nrow=6,ncol=8)
maxBen <- c(100, 2, 1, 100, 100, 100)
maxHarm <- c(100, 100, 100, 2, 1, 0)

for (i in 1:6){
  result.noBV[i,] <- boundsNoCov_res(ordinalScale, YTd, YCd, maxBen[i], maxHarm[i])
  result.BV[i,] <- boundsCov_res(ordinalScale, YT1d, YC1d, YT2d, YC2d, maxBen[i], maxHarm[i])
}

result.pop <- cbind(result.noBV,result.BV[,1:2])
result.subpop <- result.BV[,3:8]

##The tables below are presented in the Supplementary Materials
result.pop <- data.frame(result.pop)
names(result.pop) <- c("l", "u", "epsilon", "l (using BV)", "u (using BV)")
row.names(result.pop) <- c("No assumptions", "Benefit at most 2", "at most 1", 
                           "Harm at most 2", " at most 1", "No Harm")
tab <- xtable(result.pop, align = "rccccc")
digits(tab) <- 2
print(tab)

result.subpop <- data.frame(result.subpop)
names(result.subpop) <- c("l", "u", "epsilon", "l", "u", "epsilon")
row.names(result.subpop) <- c("No assumptions", "Benefit at most 2", "at most 1", 
                              "Harm at most 2", " at most 1", "No Harm")
tab <- xtable(result.subpop, align = "rcccccc")
digits(tab) <- 2
print(tab)
##########################################################################################################
##Discretize RICV using a bin length of 20 mL and get bound estimates##
##########################################################################################################
levelfun <- function(x){ 
  if(x < 0){
    y = 1
  } else if(x >= 0 & x < 20){
    y = 2
  } else {
    y = 3
  } 
  return(y)
}

YCd <- sapply(YC, levelfun)
YTd <- sapply(YT, levelfun)
YC1d <- sapply(YC1, levelfun)
YT1d <- sapply(YT1, levelfun)
YC2d <- sapply(YC2, levelfun)
YT2d <- sapply(YT2, levelfun)

ordinalScale <- 1:3

result.noBV <- matrix(0,nrow=4,ncol=3)
result.BV <- matrix(NA,nrow=4,ncol=8)
maxBen <- c(100, 1, 100, 100)
maxHarm <- c(100, 100, 1, 0)

for (i in 1:4){
  result.noBV[i,] <- boundsNoCov_res(ordinalScale, YTd, YCd, maxBen[i], maxHarm[i])
  result.BV[i,] <- boundsCov_res(ordinalScale, YT1d, YC1d, YT2d, YC2d, maxBen[i], maxHarm[i])
}

result.pop <- cbind(result.noBV,result.BV[,1:2])
result.subpop <- result.BV[,3:8]

##The tables below are presented in the Supplementary Materials
result.pop <- data.frame(result.pop)
names(result.pop) <- c("l", "u", "epsilon", "l (using BV)", "u (using BV)")
row.names(result.pop) <- c("No assumptions", "Benefit at most 1", 
                           "Harm at most 1", "No Harm")
tab <- xtable(result.pop, align = "rccccc")
digits(tab) <- 2
print(tab)

result.subpop <- data.frame(result.subpop)
names(result.subpop) <- c("l", "u", "epsilon", "l", "u", "epsilon")
row.names(result.subpop) <- c("No assumptions", "Benefit at most 1", 
                              "Harm at most 1", "No Harm")
tab <- xtable(result.subpop, align = "rcccccc")
digits(tab) <- 2
print(tab)

