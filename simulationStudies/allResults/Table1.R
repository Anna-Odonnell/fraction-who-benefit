library(xtable)
SE.fun <- function(x){sqrt(mean((x-mean(x))^2))}

###################################################################################################
##RICV5##
###################################################################################################
setwd("RICV5")

trueLB <- 0.8227513
trueUB <- 0.9629630

mat1 <- matrix(NA, nrow = 3, ncol = 4)
mat2 <- matrix(NA, nrow = 3, ncol = 8)

load("n100_nbootstrap.Rdata")
mat1[1,1] <- mean(bounds[,1])-trueLB
mat1[1,2] <- mean(bounds[,2])-trueUB
mat1[1,3] <- SE.fun(bounds[,1])
mat1[1,4] <- SE.fun(bounds[,2])

mat2[1,1] <- mean(trueLB <= LB.CI_perc[,2] & trueLB >= LB.CI_perc[,1])
mat2[1,2] <- mean(trueUB <= UB.CI_perc[,2] & trueUB >= UB.CI_perc[,1])
mat2[1,5] <- mean(LB.CI_perc[,2]-LB.CI_perc[,1])
mat2[1,6] <- mean(UB.CI_perc[,2]-UB.CI_perc[,1])

load("n500_nbootstrap.Rdata")
mat1[2,1] <- mean(bounds[,1])-trueLB
mat1[2,2] <- mean(bounds[,2])-trueUB
mat1[2,3] <- SE.fun(bounds[,1])
mat1[2,4] <- SE.fun(bounds[,2])

mat2[2,1] <- mean(trueLB <= LB.CI_perc[,2] & trueLB >= LB.CI_perc[,1])
mat2[2,2] <- mean(trueUB <= UB.CI_perc[,2] & trueUB >= UB.CI_perc[,1])
mat2[2,5] <- mean(LB.CI_perc[,2]-LB.CI_perc[,1])
mat2[2,6] <- mean(UB.CI_perc[,2]-UB.CI_perc[,1])

load("n1000_nbootstrap.Rdata")
mat1[3,1] <- mean(bounds[,1])-trueLB
mat1[3,2] <- mean(bounds[,2])-trueUB
mat1[3,3] <- SE.fun(bounds[,1])
mat1[3,4] <- SE.fun(bounds[,2])

mat2[3,1] <- mean(trueLB <= LB.CI_perc[,2] & trueLB >= LB.CI_perc[,1])
mat2[3,2] <- mean(trueUB <= UB.CI_perc[,2] & trueUB >= UB.CI_perc[,1])
mat2[3,5] <- mean(LB.CI_perc[,2]-LB.CI_perc[,1])
mat2[3,6] <- mean(UB.CI_perc[,2]-UB.CI_perc[,1])

load("n100_mofn.Rdata")
mat2[1,3] <- mean(trueLB <= LB.CI_mn[,2] & trueLB >= LB.CI_mn[,1])
mat2[1,4] <- mean(trueUB <= UB.CI_mn[,2] & trueUB >= UB.CI_mn[,1])
mat2[1,7] <- mean(LB.CI_mn[,2]-LB.CI_mn[,1])
mat2[1,8] <- mean(UB.CI_mn[,2]-UB.CI_mn[,1])

load("n500_mofn.Rdata")
mat2[2,3] <- mean(trueLB <= LB.CI_mn[,2] & trueLB >= LB.CI_mn[,1])
mat2[2,4] <- mean(trueUB <= UB.CI_mn[,2] & trueUB >= UB.CI_mn[,1])
mat2[2,7] <- mean(LB.CI_mn[,2]-LB.CI_mn[,1])
mat2[2,8] <- mean(UB.CI_mn[,2]-UB.CI_mn[,1])

load("n1000_mofn.Rdata")
mat2[3,3] <- mean(trueLB <= LB.CI_mn[,2] & trueLB >= LB.CI_mn[,1])
mat2[3,4] <- mean(trueUB <= UB.CI_mn[,2] & trueUB >= UB.CI_mn[,1])
mat2[3,7] <- mean(LB.CI_mn[,2]-LB.CI_mn[,1])
mat2[3,8] <- mean(UB.CI_mn[,2]-UB.CI_mn[,1])

mat1_RICV <- mat1
mat2_RICV <- mat2

rm(mat1,mat2)

###################################################################################################
##Binary (no restrictions)##
###################################################################################################
setwd("binary_noRestrictions")

trueLB <- 0
trueUB <- 0.5

mat1 <- matrix(NA, nrow = 3, ncol = 4)
mat2 <- matrix(NA, nrow = 3, ncol = 8)

load("n100_nbootstrap.Rdata")
mat1[1,1] <- mean(bounds[,1])-trueLB
mat1[1,2] <- mean(bounds[,2])-trueUB
mat1[1,3] <- SE.fun(bounds[,1])
mat1[1,4] <- SE.fun(bounds[,2])

mat2[1,1] <- mean(trueLB <= LB.CI_perc[,2] & trueLB >= LB.CI_perc[,1])
mat2[1,2] <- mean(trueUB <= UB.CI_perc[,2] & trueUB >= UB.CI_perc[,1])
mat2[1,5] <- mean(LB.CI_perc[,2]-LB.CI_perc[,1])
mat2[1,6] <- mean(UB.CI_perc[,2]-UB.CI_perc[,1])

load("n500_nbootstrap.Rdata")
mat1[2,1] <- mean(bounds[,1])-trueLB
mat1[2,2] <- mean(bounds[,2])-trueUB
mat1[2,3] <- SE.fun(bounds[,1])
mat1[2,4] <- SE.fun(bounds[,2])

mat2[2,1] <- mean(trueLB <= LB.CI_perc[,2] & trueLB >= LB.CI_perc[,1])
mat2[2,2] <- mean(trueUB <= UB.CI_perc[,2] & trueUB >= UB.CI_perc[,1])
mat2[2,5] <- mean(LB.CI_perc[,2]-LB.CI_perc[,1])
mat2[2,6] <- mean(UB.CI_perc[,2]-UB.CI_perc[,1])

load("n1000_nbootstrap.Rdata")
mat1[3,1] <- mean(bounds[,1])-trueLB
mat1[3,2] <- mean(bounds[,2])-trueUB
mat1[3,3] <- SE.fun(bounds[,1])
mat1[3,4] <- SE.fun(bounds[,2])

mat2[3,1] <- mean(trueLB <= LB.CI_perc[,2] & trueLB >= LB.CI_perc[,1])
mat2[3,2] <- mean(trueUB <= UB.CI_perc[,2] & trueUB >= UB.CI_perc[,1])
mat2[3,5] <- mean(LB.CI_perc[,2]-LB.CI_perc[,1])
mat2[3,6] <- mean(UB.CI_perc[,2]-UB.CI_perc[,1])

load("n100_mofn.Rdata")
mat2[1,3] <- mean(trueLB <= LB.CI_mn[,2] & trueLB >= LB.CI_mn[,1])
mat2[1,4] <- mean(trueUB <= UB.CI_mn[,2] & trueUB >= UB.CI_mn[,1])
mat2[1,7] <- mean(LB.CI_mn[,2]-LB.CI_mn[,1])
mat2[1,8] <- mean(UB.CI_mn[,2]-UB.CI_mn[,1])

load("n500_mofn.Rdata")
mat2[2,3] <- mean(trueLB <= LB.CI_mn[,2] & trueLB >= LB.CI_mn[,1])
mat2[2,4] <- mean(trueUB <= UB.CI_mn[,2] & trueUB >= UB.CI_mn[,1])
mat2[2,7] <- mean(LB.CI_mn[,2]-LB.CI_mn[,1])
mat2[2,8] <- mean(UB.CI_mn[,2]-UB.CI_mn[,1])

load("n1000_mofn.Rdata")
mat2[3,3] <- mean(trueLB <= LB.CI_mn[,2] & trueLB >= LB.CI_mn[,1])
mat2[3,4] <- mean(trueUB <= UB.CI_mn[,2] & trueUB >= UB.CI_mn[,1])
mat2[3,7] <- mean(LB.CI_mn[,2]-LB.CI_mn[,1])
mat2[3,8] <- mean(UB.CI_mn[,2]-UB.CI_mn[,1])

mat1_binNoRes <- mat1
mat2_binNoRes <- mat2

rm(mat1,mat2)

###################################################################################################
##Binary (no harm)##
###################################################################################################
setwd("binary_noHarm")

trueLB <- 0
trueUB <- 0

mat1 <- matrix(NA, nrow = 3, ncol = 4)
mat2 <- matrix(NA, nrow = 3, ncol = 8)

load("n100_nbootstrap.Rdata")
mat1[1,1] <- mean(bounds[,1])-trueLB
mat1[1,2] <- mean(bounds[,2])-trueUB
mat1[1,3] <- SE.fun(bounds[,1])
mat1[1,4] <- SE.fun(bounds[,2])

mat2[1,1] <- mean(trueLB <= LB.CI_perc[,2] & trueLB >= LB.CI_perc[,1])
mat2[1,2] <- mean(trueUB <= UB.CI_perc[,2] & trueUB >= UB.CI_perc[,1])
mat2[1,5] <- mean(LB.CI_perc[,2]-LB.CI_perc[,1])
mat2[1,6] <- mean(UB.CI_perc[,2]-UB.CI_perc[,1])

load("n500_nbootstrap.Rdata")
mat1[2,1] <- mean(bounds[,1])-trueLB
mat1[2,2] <- mean(bounds[,2])-trueUB
mat1[2,3] <- SE.fun(bounds[,1])
mat1[2,4] <- SE.fun(bounds[,2])

mat2[2,1] <- mean(trueLB <= LB.CI_perc[,2] & trueLB >= LB.CI_perc[,1])
mat2[2,2] <- mean(trueUB <= UB.CI_perc[,2] & trueUB >= UB.CI_perc[,1])
mat2[2,5] <- mean(LB.CI_perc[,2]-LB.CI_perc[,1])
mat2[2,6] <- mean(UB.CI_perc[,2]-UB.CI_perc[,1])

load("n1000_nbootstrap.Rdata")
mat1[3,1] <- mean(bounds[,1])-trueLB
mat1[3,2] <- mean(bounds[,2])-trueUB
mat1[3,3] <- SE.fun(bounds[,1])
mat1[3,4] <- SE.fun(bounds[,2])

mat2[3,1] <- mean(trueLB <= LB.CI_perc[,2] & trueLB >= LB.CI_perc[,1])
mat2[3,2] <- mean(trueUB <= UB.CI_perc[,2] & trueUB >= UB.CI_perc[,1])
mat2[3,5] <- mean(LB.CI_perc[,2]-LB.CI_perc[,1])
mat2[3,6] <- mean(UB.CI_perc[,2]-UB.CI_perc[,1])

load("n100_mofn.Rdata")
mat2[1,3] <- mean(trueLB <= LB.CI_mn[,2] & trueLB >= LB.CI_mn[,1])
mat2[1,4] <- mean(trueUB <= UB.CI_mn[,2] & trueUB >= UB.CI_mn[,1])
mat2[1,7] <- mean(LB.CI_mn[,2]-LB.CI_mn[,1])
mat2[1,8] <- mean(UB.CI_mn[,2]-UB.CI_mn[,1])

load("n500_mofn.Rdata")
mat2[2,3] <- mean(trueLB <= LB.CI_mn[,2] & trueLB >= LB.CI_mn[,1])
mat2[2,4] <- mean(trueUB <= UB.CI_mn[,2] & trueUB >= UB.CI_mn[,1])
mat2[2,7] <- mean(LB.CI_mn[,2]-LB.CI_mn[,1])
mat2[2,8] <- mean(UB.CI_mn[,2]-UB.CI_mn[,1])

load("n1000_mofn.Rdata")
mat2[3,3] <- mean(trueLB <= LB.CI_mn[,2] & trueLB >= LB.CI_mn[,1])
mat2[3,4] <- mean(trueUB <= UB.CI_mn[,2] & trueUB >= UB.CI_mn[,1])
mat2[3,7] <- mean(LB.CI_mn[,2]-LB.CI_mn[,1])
mat2[3,8] <- mean(UB.CI_mn[,2]-UB.CI_mn[,1])

mat1_binNoHarm <- mat1
mat2_binNoHarm <- mat2

rm(mat1,mat2)


###################################################################################################
##Concatenate the results##
###################################################################################################
mat1 <- rbind(mat1_RICV,mat1_binNoRes,mat1_binNoHarm)
mat2 <- rbind(mat2_RICV,mat2_binNoRes,mat2_binNoHarm)
mat1 <- data.frame(mat1)
mat2 <- data.frame(mat2)
table1A <- xtable(mat1, align = "ccccc")
table1B <- xtable(mat2, align = "ccccccccc")
digits(table1A) <- 3
digits(table1B) <- 3
print(table1A)
print(table1B)
