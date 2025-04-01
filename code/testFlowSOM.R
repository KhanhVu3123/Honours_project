require(FlowSOM)
require(flowCore)
require(dplyr)

library(mvtnorm)

sigma_6D <- rep(1,36)
i <-1 
while(i <= 36){
  sigma_6D[i] <- 2
  i <- i + 7
}
sigma_6D <- matrix(sigma_6D, nrow = 6, ncol = 6)

dat6D_1 <- rmvnorm(100000, mean = rep(0,6), sigma_6D)
dat6D_2 <- rmvnorm(100000, mean = c(10, rep(0,5)), sigma_6D)
dat6D_3 <- rmvnorm(100000, mean = c(0,10, rep(0,4)), sigma_6D)

dat6D <- rbind(dat6D_1, dat6D_2,dat6D_3)
colnames(dat6D) <- c("V1", "V2","V3","V4","V5","V6")

colnames(dat6D_1) <- c("V1", "V2","V3","V4","V5","V6")

fSOM_test <- FlowSOM::ReadInput(as.matrix(dat6D_1), transform = F, scale = F)
fSOM_test <- FlowSOM::BuildSOM(fSOM_test, colsToUse = c(1,2,3,4,5,6))
fSOM_test <- FlowSOM::BuildMST(fSOM_test)

FlowSOM::PlotStars(fSOM_test)

