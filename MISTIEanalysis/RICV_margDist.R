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
##Discretize RICV using a bin length of 5 mL and get marginal distribution plots##
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

nC <- length(YCd)
nT <- length(YTd)
nC1 <- length(YC1d)
nT1 <- length(YT1d)
nC2 <- length(YC2d)
nT2 <- length(YT2d)

##Counting sample sizes in control and treatment groups 
nLevels <- length(ordinalScale)
countC <- rep(0,nLevels)
countT <- rep(0,nLevels)
countC1 <- rep(0,nLevels)
countT1 <- rep(0,nLevels)
countC2 <- rep(0,nLevels)
countT2 <- rep(0,nLevels)

for(i in 1:nLevels){
  countC[i] <- sum(YCd==i)
  countT[i] <- sum(YTd==i)
  countC1[i] <- sum(YC1d==i)
  countT1[i] <- sum(YT1d==i)
  countC2[i] <- sum(YC2d==i)
  countT2[i] <- sum(YT2d==i)
}

pdf("RICV5_margDist.pdf")
par(mfcol=c(3,1))
mat <- rbind(countT, countC)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Total Population", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV5", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 3, y = 0.9, labels = bquote(n[T]~"="~.(nT)~","~n[C]~"="~.(nC)),cex=2)
legend(x=12,y=1,legend = rownames(mat), fill = c("red", "blue"),cex = 1.5)

mat <- rbind(countT1, countC1)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with low baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV5", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 3, y = 0.9, labels = bquote(n[T]~"="~.(nT1)~","~n[C]~"="~.(nC1)),cex=2)

mat <- rbind(countT2, countC2)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with high baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV5", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 3, y = 0.9, labels = bquote(n[T]~"="~.(nT2)~","~n[C]~"="~.(nC2)),cex=2)
dev.off()
rm(YCd,YTd,YC1d,YC2d,YT1d,YT2d)
##########################################################################################################
##Discretize RICV using a bin length of 2 mL and get marginal distribution plots##
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

nC <- length(YCd)
nT <- length(YTd)
nC1 <- length(YC1d)
nT1 <- length(YT1d)
nC2 <- length(YC2d)
nT2 <- length(YT2d)

##Counting sample sizes in control and treatment groups 
nLevels <- length(ordinalScale)
countC <- rep(0,nLevels)
countT <- rep(0,nLevels)
countC1 <- rep(0,nLevels)
countT1 <- rep(0,nLevels)
countC2 <- rep(0,nLevels)
countT2 <- rep(0,nLevels)

for(i in 1:nLevels){
  countC[i] <- sum(YCd==i)
  countT[i] <- sum(YTd==i)
  countC1[i] <- sum(YC1d==i)
  countT1[i] <- sum(YT1d==i)
  countC2[i] <- sum(YC2d==i)
  countT2[i] <- sum(YT2d==i)
}

pdf("RICV2_margDist.pdf")
par(mfcol=c(3,1))
mat <- rbind(countT, countC)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Total Population", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV2", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 5, y = 0.9, labels = bquote(n[T]~"="~.(nT)~","~n[C]~"="~.(nC)),cex=2)
legend(x=12,y=1,legend = rownames(mat), fill = c("red", "blue"),cex = 1.5)

mat <- rbind(countT1, countC1)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with low baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV2", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 5, y = 0.9, labels = bquote(n[T]~"="~.(nT1)~","~n[C]~"="~.(nC1)),cex=2)

mat <- rbind(countT2, countC2)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with high baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV2", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 5, y = 0.9, labels = bquote(n[T]~"="~.(nT2)~","~n[C]~"="~.(nC2)),cex=2)

dev.off()

rm(YCd,YTd,YC1d,YC2d,YT1d,YT2d)
##########################################################################################################
##Discretize RICV using a bin length of 10 mL and get marginal distribution plots##
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

nC <- length(YCd)
nT <- length(YTd)
nC1 <- length(YC1d)
nT1 <- length(YT1d)
nC2 <- length(YC2d)
nT2 <- length(YT2d)

##Counting sample sizes in control and treatment groups 
nLevels <- length(ordinalScale)
countC <- rep(0,nLevels)
countT <- rep(0,nLevels)
countC1 <- rep(0,nLevels)
countT1 <- rep(0,nLevels)
countC2 <- rep(0,nLevels)
countT2 <- rep(0,nLevels)

for(i in 1:nLevels){
  countC[i] <- sum(YCd==i)
  countT[i] <- sum(YTd==i)
  countC1[i] <- sum(YC1d==i)
  countT1[i] <- sum(YT1d==i)
  countC2[i] <- sum(YC2d==i)
  countT2[i] <- sum(YT2d==i)
}

pdf("RICV10_margDist.pdf")
par(mfcol=c(3,1))
mat <- rbind(countT, countC)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Total Population", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV10", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 2.3, y = 0.9, labels = bquote(n[T]~"="~.(nT)~","~n[C]~"="~.(nC)),cex=2)
legend(x=7,y=1,legend = rownames(mat), fill = c("red", "blue"),cex = 1.5)

mat <- rbind(countT1, countC1)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with low baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV10", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 2.3, y = 0.9, labels = bquote(n[T]~"="~.(nT1)~","~n[C]~"="~.(nC1)),cex=2)

mat <- rbind(countT2, countC2)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with high baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV10", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 2.3, y = 0.9, labels = bquote(n[T]~"="~.(nT2)~","~n[C]~"="~.(nC2)),cex=2)
dev.off()

rm(YCd,YTd,YC1d,YC2d,YT1d,YT2d)
##########################################################################################################
##Discretize RICV using a bin length of 20 mL and get marginal distribution plots##
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

nC <- length(YCd)
nT <- length(YTd)
nC1 <- length(YC1d)
nT1 <- length(YT1d)
nC2 <- length(YC2d)
nT2 <- length(YT2d)

##Counting sample sizes in control and treatment groups 
nLevels <- length(ordinalScale)
countC <- rep(0,nLevels)
countT <- rep(0,nLevels)
countC1 <- rep(0,nLevels)
countT1 <- rep(0,nLevels)
countC2 <- rep(0,nLevels)
countT2 <- rep(0,nLevels)

for(i in 1:nLevels){
  countC[i] <- sum(YCd==i)
  countT[i] <- sum(YTd==i)
  countC1[i] <- sum(YC1d==i)
  countT1[i] <- sum(YT1d==i)
  countC2[i] <- sum(YC2d==i)
  countT2[i] <- sum(YT2d==i)
}

pdf("RICV20_margDist.pdf")
par(mfcol=c(3,1))
mat <- rbind(countT, countC)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Total Population", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV20", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 1.95, y = 0.9, labels = bquote(n[T]~"="~.(nT)~","~n[C]~"="~.(nC)),cex=2)
legend(x=3.2,y=1,legend = rownames(mat), fill = c("red", "blue"),cex = 1.5)

mat <- rbind(countT1, countC1)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with low baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV20", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 1.95, y = 0.9, labels = bquote(n[T]~"="~.(nT1)~","~n[C]~"="~.(nC1)),cex=2)

mat <- rbind(countT2, countC2)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with high baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV20", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 1.95, y = 0.9, labels = bquote(n[T]~"="~.(nT2)~","~n[C]~"="~.(nC2)),cex=2)

dev.off()

rm(YCd,YTd,YC1d,YC2d,YT1d,YT2d)
