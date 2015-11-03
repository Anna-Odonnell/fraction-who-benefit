rm(list=ls())


##The files below were created by changeInVol.R
YC <- unlist(read.table("YCd.txt"))
YT <- unlist(read.table("YTd.txt"))
YC1 <- unlist(read.table("YC1d.txt"))
YT1 <- unlist(read.table("YT1d.txt"))
YC2 <- unlist(read.table("YC2d.txt"))
YT2 <- unlist(read.table("YT2d.txt"))

##Counting sample sizes in control and treatment groups 
nLevels <- 6
countC <- rep(0,nLevels)
countT <- rep(0,nLevels)
countC1 <- rep(0,nLevels)
countT1 <- rep(0,nLevels)
countC2 <- rep(0,nLevels)
countT2 <- rep(0,nLevels)

for(i in 1:nLevels){
  countC[i] <- sum(YC==i)
  countT[i] <- sum(YT==i)
  countC1[i] <- sum(YC1==i)
  countT1[i] <- sum(YT1==i)
  countC2[i] <- sum(YC2==i)
  countT2[i] <- sum(YT2==i)
}




##########################################################################################################
##Make Web Figure 1##
##########################################################################################################
nC <- length(YC)
nT <- length(YT)
nC1 <- length(YC1)
nT1 <- length(YT1)
nC2 <- length(YC2)
nT2 <- length(YT2)

pdf("Figure1Supp.pdf")
par(mfcol=c(3,1))
mat <- rbind(countT, countC)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Total Population", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 2.8, y = 0.9, labels = bquote(n[T]~"="~.(nT)~","~n[C]~"="~.(nC)),cex=2)
legend(x=12,y=1,legend = rownames(mat), fill = c("red", "blue"),cex = 1.5)

mat <- rbind(countT1, countC1)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with low baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 2.8, y = 0.9, labels = bquote(n[T]~"="~.(nT1)~","~n[C]~"="~.(nC1)),cex=2)

mat <- rbind(countT2, countC2)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with high baseline clot volume", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 2.8, y = 0.9, labels = bquote(n[T]~"="~.(nT2)~","~n[C]~"="~.(nC2)),cex=2)
dev.off()

