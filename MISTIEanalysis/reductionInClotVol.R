##IN THIS CODE:
##create Figure 5
##estimate bounds for reduction in clot volume

##inclusion criteria:
##randomized, have reduction in clot volume, have baseline clot volume
##treatment and control patients from MISTIE II
##control patients from ICES group

##########################################################################################################
##Obtain YC, YT, YC1, YC2, YT1, YT2##
##########################################################################################################


rm(list=ls())

##Load the Rdata file with the MISTIE II dataset
##It includes three data frames: baselineData, mRSData, and clotData
##Each of these data frames has 96 people (42 medical, 54 surgical)
##All inclusion criteria have been met except have reduction in clot volume 
##and have baseline clot volume

##We can ignore mrsData
rm(mrsData)

##Merge baselineData and clotData
data <- merge(baselineData, clotData, by = "patientName")

rm(baselineData, clotData)

##Check that merge was successful
all(data$Group_Assigned == data$group_assigned) ##true

##Everyone has reduction in clot volume
sum(is.na(data$Pre_Rand_ICH_Volume_RC)) ##0
sum(is.na(data$eot_ich_9_13)) ##0

data$volchange <- data$Pre_Rand_ICH_Volume_RC - data$eot_ich_9_13

##There will be 96 patients in this analysis.

##Focus on the variables of interest
data <- data.frame(data$patientName, data$Group_Assigned, data$volchange, data$Pre_Rand_ICH_Volume_RC)
names(data) <- c("id", "group", "volChange", "baseVol")

control <- subset(data, group == "Medical") ##42
tmt <- subset(data, group == "Surgical") ##54

YT <- tmt$volChange
YC <- control$volChange


m <- median(data$baseVol) ##43.204999925

control1 <- subset(control, baseVol < m) ##22
tmt1 <- subset(tmt, baseVol < m) ##26
control2 <- subset(control, baseVol >= m) ##20
tmt2 <- subset(tmt, baseVol >= m) ##28

YT1 <- tmt1$volChange
YT2 <- tmt2$volChange
YC1 <- control1$volChange
YC2 <- control2$volChange


##################################################################################
##Code for Figure 5##
##################################################################################

##PLOT OF EMPIRICAL CUMULATIVE DISTRIBUTION FUNCTIONS
YT_sort <- c(-20, sort(unique(YT)),80)
YC_sort <- c(-20, sort(unique(YC)),80)
lT <- length(YT_sort)
lC <- length(YC_sort)
FT <- rep(NA, lT)
FC <- rep(NA, lC)
for(i in 1:lT){
  FT[i] <- sum(YT <= YT_sort[i])
}
for(j in 1:lC){
  FC[j] <- sum(YC <= YC_sort[j])
}
FT <- FT/length(YT)
FC <- FC/length(YC)


YT1_sort <- c(-20, sort(unique(YT1)),80)
YC1_sort <- c(-20, sort(unique(YC1)),80)
lT1 <- length(YT1_sort)
lC1 <- length(YC1_sort)
FT1 <- rep(NA, lT1)
FC1 <- rep(NA, lC1)
for(i in 1:lT1){
  FT1[i] <- sum(YT1 <= YT1_sort[i])
}
for(j in 1:lC1){
  FC1[j] <- sum(YC1 <= YC1_sort[j])
}
FT1 <- FT1/length(YT1)
FC1 <- FC1/length(YC1)

YT2_sort <- c(-20, sort(unique(YT2)),80)
YC2_sort <- c(-20, sort(unique(YC2)),80)
lT2 <- length(YT2_sort)
lC2 <- length(YC2_sort)
FT2 <- rep(NA, lT2)
FC2 <- rep(NA, lC2)
for(i in 1:lT2){
  FT2[i] <- sum(YT2 <= YT2_sort[i])
}
for(j in 1:lC2){
  FC2[j] <- sum(YC2 <= YC2_sort[j])
}
FT2 <- FT2/length(YT2)
FC2 <- FC2/length(YC2)


pdf("Figure5.pdf")
par(mfcol=c(3,1))
plot(YT_sort, FT, type = "s", xlim = c(-20,80), ylim = c(0,1), xlab = "Reduction in clot volume (mL)", ylab = NA, col = "red", main = "Total population")
par(new = TRUE)
plot(YC_sort, FC, type = "s", axes = F, xlim = c(-20,80), ylim = c(0,1), xlab = NA, ylab = NA, col = "blue")
plot(YT1_sort, FT1, type = "s", xlim = c(-20,80), ylim = c(0,1), xlab = "Reduction in clot volume (mL)", ylab = NA, col = "red", main = "Subpopulation with low baseline clot volume")
par(new = TRUE)
plot(YC1_sort, FC1, type = "s", axes = F, xlim = c(-20,80), ylim = c(0,1), xlab = NA, ylab = NA, col = "blue")
plot(YT2_sort, FT2, type = "s", xlim = c(-20,80), ylim = c(0,1), xlab = "Reduction in clot volume (mL)", ylab = NA, col = "red", main = "Subpopulation with high baseline clot volume")
par(new = TRUE)
plot(YC2_sort, FC2, type = "s", axes = F, xlim = c(-20,80), ylim = c(0,1), xlab = NA, ylab = NA, col = "blue")
dev.off()


##########################################################################################################
##Calculate bound estimates##
##########################################################################################################


write.table(YC, file = "YC.txt", row.names = FALSE, col.names = FALSE)
write.table(YT, file = "YT.txt", row.names = FALSE, col.names = FALSE)
write.table(YC1, file = "YC1.txt", row.names = FALSE, col.names = FALSE)
write.table(YT1, file = "YT1.txt", row.names = FALSE, col.names = FALSE)
write.table(YC2, file = "YC2.txt", row.names = FALSE, col.names = FALSE)
write.table(YT2, file = "YT2.txt", row.names = FALSE, col.names = FALSE)

##Open MATLAB and run the rest of this code in MATLAB. This code uses the bound
##functions, which are also provided on the repository. Make sure your working directory
##is the one in which you have saved these functions and "YC.txt", "YC1.txt", etc.


YT = importdata('YT.txt');
YC = importdata('YC.txt');
YT1 = importdata('YT1.txt');
YT2 = importdata('YT2.txt');
YC1 = importdata('YC1.txt');
YC2 = importdata('YC2.txt');



[l,u,lflag,uflag] = boundsNoCov(YT, YC)
[l,u,lflag, uflag] = boundsNoCov_resHarm(YT, YC, 0)
[l,u,lflag, uflag] = boundsNoCov_resHarm(YT, YC, 1)
[l,u,lflag, uflag] = boundsNoCov_resHarm(YT, YC, 2)
[l,u,lflag, uflag] = boundsNoCov_resHarm(YT, YC, 3)
[l,u,lflag, uflag] = boundsNoCov_resHarm(YT, YC, 4)
[l,u,lflag, uflag] = boundsNoCov_resHarm(YT, YC, 5)

[l, u, l1, u1, l2, u2, flag] = boundsCov(YT1, YC1, YT2, YC2)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 0, 0)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 1, 1)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 2, 2)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 3, 3)
[l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, 4, 4)




