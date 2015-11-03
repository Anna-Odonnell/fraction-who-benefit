##The analysis below includes treatment and control patients from MISTIE II, and control patients from ICES group.
##Patients are included only if they had random treatment assignment, non-missing 180-day mRS score, 
##and non-missing baseline NIHSS score

##########################################################################################################
##Gather the data for 180-day mRS score##
##########################################################################################################
rm(list=ls())

load("MISTIEIIdata.Rdata")
##Each of these datasets (i.e., mRSData, clotData, baselineData) has 96 people (42 medical, 54 surgical)
##All inclusion criteria have been met except non-missing 180-day mRS score and baseline NIHSS score

##We only need mRSData and baselineData right now
rm(clotData)

##Also, we want to focus on the data at 180-days
mrs180 <- subset(mrsData, Follow_up_Visit == 180) 
rm(mrsData)
length(unique(mrs180$patientName)) ##96

##Merge mrs180 with baselineData
data <- merge(baselineData, mrs180, by = "patientName")

##Check that merge was successful
all(data$Group_Assigned.x == data$Group_Assigned.y) ##true

##Drop redundant columns
drops <- c("RunIn_Randomized.y","Group_Assigned.y")
data <- data[,!(names(data) %in% drops)]

rm(drops, mrs180, baselineData)

##Only include people with baseline NIHSS and 180-day mRS score
data$patientName[is.na(data$Enrollment_NIHSS_Total)] 
data <- subset(data, !is.na(Enrollment_NIHSS_Total))
nrow(data) ##95

sum(is.na(data$rankin_score)) ##6 people missing rankin_score
data$patientName[is.na(data$rankin_score)] 
data <- subset(data, !is.na(rankin_score)) 
nrow(data) ##89

##There will be 89 people total in this analysis

control <- subset(data, Group_Assigned.x == "Medical") ##37
treatment <- subset(data, Group_Assigned.x == "Surgical") ##52

YC <- control$rankin_score 
YT <- treatment$rankin_score

control1 <- subset(control, Enrollment_NIHSS_Total <= 20) ##17
control2 <- subset(control, Enrollment_NIHSS_Total > 20) ##20
tmt1 <- subset(treatment, Enrollment_NIHSS_Total <= 20) ##21
tmt2 <- subset(treatment, Enrollment_NIHSS_Total > 20) ##31

YC1 <- control1$rankin_score
YC2 <- control2$rankin_score
YT1 <- tmt1$rankin_score
YT2 <- tmt2$rankin_score

rm(control1, control2, tmt1, tmt2, treatment, control, data)


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

YC <- sapply(YC, levelfun)
YT <- sapply(YT, levelfun)
YC1 <- sapply(YC1, levelfun)
YT1 <- sapply(YT1, levelfun)
YC2 <- sapply(YC2, levelfun)
YT2 <- sapply(YT2, levelfun)


write.table(YC, file = "YC.txt", col.names = FALSE, row.names = FALSE)
write.table(YT, file = "YT.txt", col.names = FALSE, row.names = FALSE)
write.table(YC1, file = "YC1.txt", col.names = FALSE, row.names = FALSE)
write.table(YC2, file = "YC2.txt", col.names = FALSE, row.names = FALSE)
write.table(YT1, file = "YT1.txt", col.names = FALSE, row.names = FALSE)
write.table(YT2, file = "YT2.txt", col.names = FALSE, row.names = FALSE)

##########################################################################################################
##Compute the bound estimates (we use MATLAB for this step)##
##########################################################################################################
##In MATLAB


YT = importdata('YT.txt');
YC = importdata('YC.txt');
YT1 = importdata('YT1.txt');
YT2 = importdata('YT2.txt');
YC1 = importdata('YC1.txt');
YC2 = importdata('YC2.txt');


res1 = zeros(12,5);
res2 = zeros(12,6);

[~,~,l,u,eps] = boundsNoCov_res(YT, YC, 100, 100);
res1(1,1:3)=[l,u,eps];
[l, u, l1, u1, eps1, l2, u2, eps2] = boundsCov_res(YT1, YT2, YC1, YC2, 100, 100);
res1(1,4:5)=[l,u];
res2(1,:)=[l1,u1,eps1,l2,u2,eps2];

for i = 0:5
[~,~,l,u,eps] = boundsNoCov_res(YT, YC, 100, i);
res1(end-i,1:3)=[l,u,eps];
[l, u, l1, u1, eps1, l2, u2, eps2] = boundsCov_res(YT1, YT2, YC1, YC2, 100, i);
res1(end-i,4:5)=[l,u];
res2(end-i,:)=[l1,u1,eps1,l2,u2,eps2];
end

for i = 1:5
[~,~,l,u,eps] = boundsNoCov_res(YT, YC, i, 100);
res1(end-5-i,1:3)=[l,u,eps];
[l, u, l1, u1, eps1, l2, u2, eps2] = boundsCov_res(YT1, YT2, YC1, YC2, i, 100);
res1(end-5-i,4:5)=[l,u];
res2(end-5-i,:)=[l1,u1,eps1,l2,u2,eps2];
end

save('res.mat', 'res1', 'res2')

##########################################################################################################
##Make LaTeX tables of the bound estimates (Web Table 2)##
##########################################################################################################
##back in R

library(R.matlab)
library(xtable)
res <- readMat("res.mat")

res.pop <- res[[1]]
res.pop <- data.frame(res.pop)
names(res.pop) <- c("l", "u", "epsilon", "l (using BV)", "u (using BV)")
row.names(res.pop) <- c("No assumptions", "Benefit at most 5", "at most 4", "at most 3", "at most 2", "at most 1", 
                        "Harm at most 5", " at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(res.pop, align = "rccccc")
digits(tab) <- 2
print(tab)

res.subpop <- res[[2]]
res.subpop <- data.frame(res.subpop)
names(res.subpop) <- c("l", "u", "epsilon", "l", "u", "epsilon")
row.names(res.subpop) <- c("No assumptions", "Benefit at most 5", "at most 4", "at most 3", "at most 2", "at most 1", 
                           "Harm at most 5", " at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(res.subpop, align = "rcccccc")
digits(tab) <- 2
print(tab)

