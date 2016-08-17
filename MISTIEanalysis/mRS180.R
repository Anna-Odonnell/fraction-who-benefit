##The analysis below includes treatment and control patients from MISTIE II, and control patients from ICES group.
##Patients are included only if they had random treatment assignment, non-missing 180-day mRS score, 
##and non-missing baseline NIHSS score

##########################################################################################################
##Gather the data for 180-day mRS score##
##########################################################################################################
rm(list=ls())
setwd("~/Dropbox/research/data/MISTIEII")
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

##########################################################################################################
##Compute the bound estimates##
##########################################################################################################
ordinalScale <- 1:7

result.noBV <- matrix(0,nrow=12,ncol=3)
result.BV <- matrix(NA,nrow=12,ncol=8)
maxBen <- c(100, 5, 4, 3, 2, 1, 100, 100, 100, 100, 100, 100)
maxHarm <- c(100, 100, 100, 100, 100, 100, 5, 4, 3, 2, 1, 0)

for (i in 1:12){
  result.noBV[i,] <- boundsNoCov_res(ordinalScale, YT, YC, maxBen[i], maxHarm[i])
  result.BV[i,] <- boundsCov_res(ordinalScale, YT1, YC1, YT2, YC2, maxBen[i], maxHarm[i])
}

result.pop <- cbind(result.noBV,result.BV[,1:2])
result.subpop <- result.BV[,3:8]
##########################################################################################################
##Make LaTeX tables of the bound estimates (for Supplementary Materials)##
##########################################################################################################
library(xtable)

result.pop <- data.frame(result.pop)
names(result.pop) <- c("l", "u", "epsilon", "l (using BV)", "u (using BV)")
row.names(result.pop) <- c("No assumptions", "Benefit at most 5", "at most 4", "at most 3", "at most 2", "at most 1", 
                           "Harm at most 5", " at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(result.pop, align = "rccccc")
digits(tab) <- 2
print(tab)

result.subpop <- data.frame(result.subpop)
names(result.subpop) <- c("l", "u", "epsilon", "l", "u", "epsilon")
row.names(result.subpop) <- c("No assumptions", "Benefit at most 5", "at most 4", "at most 3", "at most 2", "at most 1", 
                              "Harm at most 5", " at most 4", " at most 3", " at most 2", " at most 1", "No Harm")
tab <- xtable(result.subpop, align = "rcccccc")
digits(tab) <- 2
print(tab)




##########################################################################################################
##Plot marginal distribution of 30-day mRS under treatment and control,
##with and without stratifying by BV
##########################################################################################################
nC <- length(YC)
nT <- length(YT)
nC1 <- length(YC1)
nT1 <- length(YT1)
nC2 <- length(YC2)
nT2 <- length(YT2)


countC <- rep(0,7)
countT <- rep(0,7)
countC1 <- rep(0,7)
countT1 <- rep(0,7)
countC2 <- rep(0,7)
countT2 <- rep(0,7)

for(i in 1:7){
  countC[i] <- sum(YC==i)
  countT[i] <- sum(YT==i)
  countC1[i] <- sum(YC1==i)
  countT1[i] <- sum(YT1==i)
  countC2[i] <- sum(YC2==i)
  countT2[i] <- sum(YT2==i)
}

pdf("mRS180_margDist.pdf")
par(mfcol=c(3,1))
mat <- rbind(countT, countC)
colnames(mat) <- 1:7
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Total Population", col=c("gray45","gray85"), beside=TRUE, ylim = c(0,1), yaxt="n",
        xlab = "mRS score", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 3.3, y = 0.9, labels = bquote(n[T]~"="~.(nT)~","~n[C]~"="~.(nC)),cex=2)
legend("topright",legend = rownames(mat), fill = c("gray45", "gray85"), cex = 1.75)

mat <- rbind(countT1, countC1)
colnames(mat) <- 1:7
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with non-severe stroke", col=c("gray45","gray85"), beside=TRUE, ylim = c(0,1), yaxt="n",
        xlab = "mRS score", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 3.3, y = 0.9, labels = bquote(n[T]~"="~.(nT1)~","~n[C]~"="~.(nC1)),cex=2)

mat <- rbind(countT2, countC2)
colnames(mat) <- 1:7
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with severe stroke", col=c("gray45","gray85"), beside=TRUE, ylim = c(0,1), yaxt="n",
        xlab = "mRS score", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 3.3, y = 0.9, labels = bquote(n[T]~"="~.(nT2)~","~n[C]~"="~.(nC2)),cex=2)
dev.off()



