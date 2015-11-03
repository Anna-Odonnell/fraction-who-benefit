rm(list=ls())




##The files below were created in mRS30.R
YC <- unlist(read.table("YC.txt"))
YT <- unlist(read.table("YT.txt"))
YC1 <- unlist(read.table("YC1.txt"))
YT1 <- unlist(read.table("YT1.txt"))
YC2 <- unlist(read.table("YC2.txt"))
YT2 <- unlist(read.table("YT2.txt"))


##Counting sample sizes in control and treatment groups 
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





##########################################################################################################
##Make Figure 2a##
##########################################################################################################

nC <- length(YC)
nT <- length(YT)
nC1 <- length(YC1)
nT1 <- length(YT1)
nC2 <- length(YC2)
nT2 <- length(YT2)

pdf("Figure2a.pdf")
par(mfcol=c(3,1))
mat <- rbind(countT, countC)
colnames(mat) <- 1:7
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Total Population", col=c("red","blue"), beside=TRUE, ylim = c(0,1), yaxt = "n",
        xlab = "mRS score", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 3.3, y = 0.9, labels = bquote(n[T]~"="~.(nT)~","~n[C]~"="~.(nC)),cex=2)
legend("topright",legend = rownames(mat), fill = c("red", "blue"), cex = 1.75)

mat <- rbind(countT1, countC1)
colnames(mat) <- 1:7
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with non-severe stroke", col=c("red","blue"), beside=TRUE, ylim = c(0,1), yaxt="n",
        xlab = "mRS score", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 3.3, y = 0.9, labels = bquote(n[T]~"="~.(nT1)~","~n[C]~"="~.(nC1)),cex=2)


mat <- rbind(countT2, countC2)
colnames(mat) <- 1:7
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="Subpopulation with severe stroke", col=c("red","blue"), beside=TRUE, ylim = c(0,1), yaxt="n",
        xlab = "mRS score", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","","","","","1"),las = 2, cex.axis = 2)
text(x = 3.3, y = 0.9, labels = bquote(n[T]~"="~.(nT2)~","~n[C]~"="~.(nC2)),cex=2)
dev.off()




