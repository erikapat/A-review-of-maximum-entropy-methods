#it is necessary install the function
#install.packages("Bolstad", lib="/my/own/R-packages/")
library(Bolstad)
library(pracma)

#gradient of the objetive function
#a vector with k elements with the value of the gradient
#trap is used to indicate the integration method

function(lambda, alpha, mu, correction, trap = 1)
{
 #global=0.001 #better with 0.001

 #verifications
 # lambda,alpha size bigger than 0
 if ( length(alpha) <= 0 | length(lambda) <= 0 | length(mu) <= 0){
   return(1000)
 }else if (length(alpha) != length(lambda) | length(alpha) != length(mu) | length(mu) != length(lambda)) {
   return(1001)
}
 #the objective function have two parts: an integral and escalar product
 #integrate density from 0 to 1
 estimate <- NULL
 y <- seq(from = 0, to = 1, by = global) 

 #value of Z
 #call the sum_of exponential
 #sum_e <- source("sum_exp.R")$value
 noZ   <- sum_exp(lambda, alpha, y)
 if (length(noZ) == 1) if (noZ == 1001) return(1001)

 if (trap == 0){
  Z <- sintegral(y, noZ)$value #Z 
 } else if (trap == 1){
  Z <- trapz(y, noZ) 
 }

 #exp of the sum
 sum <- sum_exp(lambda, alpha, y)
 if (length(sum) == 1 & sum==1001) {
   return(1001) # error
 }
 
 estimate <- NULL
 for (i in 1:length(alpha)){
  fy <- (-y^alpha[i])*sum
  if (trap == 0){
    int_fy <- sintegral(y,fy)$value
  } else if (trap == 1){
    int_fy <- trapz(y,fy)
  }
  #original
  if (correction == 0){
      estimate[i] <- (1/Z)*int_fy + mu[i] 
  } else if (correction == 1){
      estimate[i] <- (1/Z)*int_fy + mu[i] + 0.001*lambda[i]
  }
 }

 return(estimate)
}

#verifications
#the size of lambda, mu and alpha should be equal.