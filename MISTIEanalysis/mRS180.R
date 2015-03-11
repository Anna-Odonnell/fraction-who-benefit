##IN THIS CODE:
##create Figure 2b
##estimate bounds for 180-day mRS score

##inclusion criteria:
##randomized, have 180-day mRS score, have NIHSS score
##treatment and control patients from MISTIE II
##control patients from ICES group

##########################################################################################################
##Obtain YC, YT, YC1, YC2, YT1, YT2##
##########################################################################################################

rm(list=ls())

##Load the Rdata file with the MISTIE II dataset
##It includes three data frames: baselineData, mRSData, and clotData
##Each of these data frames has 96 people (42 medical, 54 surgical)
##All inclusion criteria have been met except have 180-day mRS score
##and have NIHSS score

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


##update scores to the scale in the paper
##Level 1 = 6; Level 2 = 5; Level 3 = 4; Level 4 = 3; Level 5 = 2; Level 6 = 1; Level 7 = 0
levelfun <- function(x){ #code mRS as levels
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

##########################################################################################################
##Code for Figure 2b##
##########################################################################################################

FT <- rep(NA, 7)
FC <- rep(NA, 7)
FT1 <- rep(NA, 7)
FC1 <- rep(NA, 7)
FT2 <- rep(NA, 7)
FC2 <- rep(NA, 7)

for(i in 1:7){
  FT[i] <- sum(YT <= i)/length(YT)
  FT1[i] <- sum(YT1 <= i)/length(YT1)
  FT2[i] <- sum(YT2 <= i)/length(YT2)
  FC[i] <- sum(YC <= i)/length(YC)
  FC1[i] <- sum(YC1 <= i)/length(YC1)
  FC2[i] <- sum(YC2 <= i)/length(YC2)
}


pdf("Figure2b.pdf")
par(mfcol=c(3,1))
plot(1:7, FC, col = "blue", type = "s", xlim = c(1,7), ylim = c(0,1), xlab = "level", ylab = "", main = "Total population")
par(new = TRUE)
plot(1:7, FT, col = "red", type = "s", xlim = c(1,7), ylim = c(0,1), xlab = "", ylab = "")

plot(1:7, FC1, col = "blue", type = "s", xlim = c(1,7), ylim = c(0,1), xlab = "level", ylab = "", main = "Subpopulation with non-severe stroke")
par(new = TRUE)
plot(1:7, FT1, col = "red", type = "s", xlim = c(1,7), ylim = c(0,1), xlab = "", ylab = "")

plot(1:7, FC2, col = "blue", type = "s", xlim = c(1,7), ylim = c(0,1), xlab = "level", ylab = "", main = "Subpopulation with severe stroke")
par(new = TRUE)
plot(1:7, FT2, col = "red", type = "s", xlim = c(1,7), ylim = c(0,1), xlab = "", ylab = "")
dev.off()

##########################################################################################################
##Calculate bound estimates##
##########################################################################################################

write.table(YC, file = "YC.txt", col.names = FALSE, row.names = FALSE)
write.table(YT, file = "YT.txt", col.names = FALSE, row.names = FALSE)
write.table(YC1, file = "YC1.txt", col.names = FALSE, row.names = FALSE)
write.table(YC2, file = "YC2.txt", col.names = FALSE, row.names = FALSE)
write.table(YT1, file = "YT1.txt", col.names = FALSE, row.names = FALSE)
write.table(YT2, file = "YT2.txt", col.names = FALSE, row.names = FALSE)

##Open MATLAB and run the rest of this code in MATLAB. This code uses the bound
##functions, which are also provided on the repository. Make sure your working directory
##is the one in which you have saved these functions and "YC.txt", "YC1.txt", etc.


YT = importdata('YT.txt');
YC = importdata('YC.txt');
YT1 = importdata('YT1.txt');
YT2 = importdata('YT2.txt');
YC1 = importdata('YC1.txt');
YC2 = importdata('YC2.txt');


[l,u,lflag,uflag] = boundsNoCov(YT,YC)
[l,u,lflag,uflag] = boundsNoCov_resBen(YT, YC, 1)
[l,u,lflag,uflag] = boundsNoCov_resBen(YT, YC, 2)
[l,u,lflag,uflag] = boundsNoCov_resBen(YT, YC, 3)
[l,u,lflag,uflag] = boundsNoCov_resBen(YT, YC, 4)
[l,u,lflag,uflag] = boundsNoCov_resBen(YT, YC, 5)

[l,u,lflag,uflag] = boundsNoCov_resHarm(YT, YC, 0)
[l,u,lflag,uflag] = boundsNoCov_resHarm(YT, YC, 1)
[l,u,lflag,uflag] = boundsNoCov_resHarm(YT, YC, 2)
[l,u,lflag,uflag] = boundsNoCov_resHarm(YT, YC, 3)
[l,u,lflag,uflag] = boundsNoCov_resHarm(YT, YC, 4)
[l,u,lflag,uflag] = boundsNoCov_resHarm(YT, YC, 5)

[l, u, l1, u1, l2, u2, flag] = boundsCov(YT1, YC1, YT2, YC2)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resBen(YT1, YT2, YC1, YC2, 1, 1)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resBen(YT1, YT2, YC1, YC2, 2, 2)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resBen(YT1, YT2, YC1, YC2, 3, 3)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resBen(YT1, YT2, YC1, YC2, 4, 4)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resBen(YT1, YT2, YC1, YC2, 5, 5)

[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 0, 0)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 1, 1)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 2, 2)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 3, 3)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 4, 4)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 5, 5)



