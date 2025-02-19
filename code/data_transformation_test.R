# This script will test 4 different data transformation method: arcsinh, biexponential (logicle), the linlog and the generalized Box-Cox transformation.
# Math will be written in the lab book, here we just implement them.

# 1. Test the different transformation method with simulated data.

# 1.1. The arcsinh method
library(ggplot2)
library(cowplot)

dat_log_norm_1 <- exp(rnorm(1000, mean = 0 , sd = 1))
hist(dat_log_norm_1)

par(mfrow = c(5,2))
for (cofactor in 1:10){
  dat_log_norm_transformed_1 <- asinh(dat_log_norm_1/ cofactor)
  hist(dat_log_norm_transformed_1, main = paste("Cofactor =", cofactor), breaks = 100)
  
}

for (cofactor in seq(10,100,10)){
  dat_log_norm_transformed_2 <- asinh(dat_log_norm_2/ cofactor)
  hist(dat_log_norm_transformed_2, main = paste("Cofactor =", cofactor), breaks = 10)
  
}


# 1.2. The linlog method

linlog <- function(y, theta){
  new_y <- c()
  for(value in y){
    if(value <= theta){
      new_value <- (value-theta)/theta + log(theta)
    }
    else{
      new_value <- log(value)
    }
    new_y <- c(new_y, new_value)
  }
  return(new_y)
}

par(mfrow = c(3,2))
hist(dat_log_norm_1)
linlog_theta_set <- c(0,1,5,10,20)

for(linlog_theta in linlog_theta_set){
  hist(linlog(dat_log_norm_1, linlog_theta))
}


dat_log_norm_2 <- exp(rnorm(1000, mean = 6 , sd = 3)) - exp(5)

par(mfrow = c(3,2))
hist(dat_log_norm_2, main = paste("Original data"))
linlog_theta_set_1 <- c(0,1,5,10,20)
linlog_theta_set_2 <- c(30,50,80,100,200)
linlog_theta_set_3 <- seq(200,300,25)

for(linlog_theta in linlog_theta_set_3){
  hist(linlog(dat_log_norm_2, linlog_theta), main = paste("Theta =", linlog_theta))
}

# 1.3. The Box-Cox transformation
boxcox_transformation <- function(y, theta){
  new_y <- c()
  for(value in y){
    new_val <- (sign(value) * abs(value)^theta - 1)/theta
    new_y <- c(new_y, new_val)
  }
  return(new_y)
}

boxcox_theta_set <- c(-5,-3,-1,1,3)
par(mfrow = c(3,2))
hist(dat_log_norm_2)
for(theta in boxcox_theta_set){
  hist(boxcox_transformation(dat_log_norm_2,theta = theta), main = paste("Theta =", theta), breaks = 50)
}

