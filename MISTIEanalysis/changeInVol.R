##Analysis below includes treatment and control patients from MISTIE II, and control patients from ICES group.
##Patients are included only if they had random treatment assignment, non-missing reduction in clot volume, 
##and non-missing baseline clot volume

rm(list=ls())

##########################################################################################################
##Gather the data about reduction in clot volume##
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


##convert continuous clot volume into a discrete variable
levelfun <- function(x){ 
  if(x < -5){
    y = 1
  } else if(x >= -5 & x < 0){
    y = 2
  } else if(x >= 0 & x < 5){
    y = 3
  } else if(x >= 5 & x < 10){
    y = 4
  } else if(x >= 10 & x < 15){
    y = 5
  } else {
    y = 6
  } 
  return(y)
}

##This is the data for RICV (marginal and stratified by low and high baseline clot volume)
YCd <- sapply(YC, levelfun)
YTd <- sapply(YT, levelfun)
YC1d <- sapply(YC1, levelfun)
YT1d <- sapply(YT1, levelfun)
YC2d <- sapply(YC2, levelfun)
YT2d <- sapply(YT2, levelfun)

##save the data
write.table(YCd, file = "YCd.txt", row.names = FALSE, col.names = FALSE)
write.table(YTd, file = "YTd.txt", row.names = FALSE, col.names = FALSE)
write.table(YC1d, file = "YC1d.txt", row.names = FALSE, col.names = FALSE)
write.table(YT1d, file = "YT1d.txt", row.names = FALSE, col.names = FALSE)
write.table(YC2d, file = "YC2d.txt", row.names = FALSE, col.names = FALSE)
write.table(YT2d, file = "YT2d.txt", row.names = FALSE, col.names = FALSE)


##########################################################################################################
##Compute the bound estimates (we use MATLAB for this step)##
##########################################################################################################

##In MATLAB

YT = importdata('YTd.txt');
YC = importdata('YCd.txt');
YT1 = importdata('YT1d.txt');
YT2 = importdata('YT2d.txt');
YC1 = importdata('YC1d.txt');
YC2 = importdata('YC2d.txt');

res1 = zeros(10,5);
res2 = zeros(10,6);

[~,~,l,u,eps] = boundsNoCov_res(YT, YC, 100, 100);
res1(1,1:3)=[l,u,eps];
[l, u, l1, u1, eps1, l2, u2, eps2] = boundsCov_res(YT1, YT2, YC1, YC2, 100, 100);
res1(1,4:5)=[l,u];
res2(1,:)=[l1,u1,eps1,l2,u2,eps2];

for i = 0:4
[~,~,l,u,eps] = boundsNoCov_res(YT, YC, 100, i);
res1(end-i,1:3)=[l,u,eps];
[l, u, l1, u1, eps1, l2, u2, eps2] = boundsCov_res(YT1, YT2, YC1, YC2, 100, i);
res1(end-i,4:5)=[l,u];
res2(end-i,:)=[l1,u1,eps1,l2,u2,eps2];
end

for i = 1:4
[~,~,l,u,eps] = boundsNoCov_res(YT, YC, i, 100);
res1(end-4-i,1:3)=[l,u,eps];
[l, u, l1, u1, eps1, l2, u2, eps2] = boundsCov_res(YT1, YT2, YC1, YC2, i, 100);
res1(end-4-i,4:5)=[l,u];
res2(end-4-i,:)=[l1,u1,eps1,l2,u2,eps2];
end

save('res.mat', 'res1', 'res2')

##########################################################################################################
##Make LaTeX tables of the bound estimates (Web Table 3)##
##########################################################################################################

##back in R

library(R.matlab)
library(xtable)

res <- readMat("res.mat")

##Web Table 4a
res.pop <- res[[1]]
res.pop <- data.frame(res.pop)
names(res.pop) <- c("l", "u", "epsilon", "l (using BV)", "u (using BV)")
row.names(res.pop) <- c("No assumptions", "Benefit at most 4", "at most 3", "at most 2", "at most 1", 
                        "Harm at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(res.pop, align = "rccccc")
digits(tab) <- 2
print(tab)

##Web Table 4b
res.subpop <- res[[2]]
res.subpop <- data.frame(res.subpop)
names(res.subpop) <- c("l", "u", "epsilon", "l", "u", "epsilon")
row.names(res.subpop) <- c("No assumptions", "Benefit at most 4", "at most 3", "at most 2", "at most 1", 
                           "Harm at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(res.subpop, align = "rcccccc")
digits(tab) <- 2
print(tab)


