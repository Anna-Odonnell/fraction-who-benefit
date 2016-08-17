##7/19/16

##The plot will have two panels

##On the left panel, show the distribution of sqrt(n)[psiHAT_l - psi_l(P_0)], under P_0, 
##for n = 100, 1000, 10000, 10^6

##On the right panel, show the distribution of sqrt(n)[psiHAT_l - psi_l(P_n)], under P_n, 
##for n = 100, 1000, 10000, 10^6

rm(list=ls())


setwd("~/Dropbox/research/github/fraction-who-benefit/nonregular")


pdf("nonregularPlot_SuppMat.pdf")
par(mfcol = c(4,2))

##Left panel
load("distributions_P0.Rdata")

h1 = hist(result100, breaks = 50, probability = TRUE, xlim = c(-3,3), ylim = c(0, .5), main = "n = 100", xlab = "")


h2 = hist(result1000, breaks = 50, probability = TRUE, xlim = c(-3,3), ylim = c(0, .5), main = "n = 1000", xlab = "")


h3 = hist(result10000, breaks = 50, probability = TRUE, xlim = c(-3,3), ylim = c(0, .5), main = "n = 10,000", xlab = "")


h4 = hist(resultMillion, breaks = 50, probability = TRUE, xlim = c(-3,3), ylim = c(0, .5), main = "n = 1,000,000", xlab = "")

curve(dnorm(x, mean=0, sd=1), col="darkblue", lwd=2, add=TRUE, yaxt="n")


##Right panel

load("distributions_Pn.Rdata")

h5 = hist(result100, breaks = 50, probability = TRUE, xlim = c(-3,3), ylim = c(0, .5), main = "n = 100", xlab = "")


h6 = hist(result1000, breaks = 50, probability = TRUE, xlim = c(-3,3), ylim = c(0, .5), main = "n = 1000", xlab = "")


h7 = hist(result10000, breaks = 50, probability = TRUE, xlim = c(-3,3), ylim = c(0, .5), main = "n = 10,000", xlab = "")


h8 = hist(resultMillion, breaks = 50, probability = TRUE, xlim = c(-3,3), ylim = c(0, .5), main = "n = 1,000,000", xlab = "")

curve(dnorm(x, mean=0, sd=1), col="darkblue", lwd=2, add=TRUE, yaxt="n")
dev.off()

##Compute the point masses for n = 10^6
load("distributions_P0.Rdata")
sum(resultMillion == 0)/length(resultMillion)

load("distributions_Pn.Rdata")
sum(resultMillion == -1)/length(resultMillion)
