#source("SME_functions.R")

Maxent_Density <- function(lambda, input, factor, minVal = 0, maxVal = 2, interVals = 30){
  
  x        <- seq(minVal/factor, maxVal/factor, length = interVals)
  densidad <- density.values(x, input, lambda)
  densidad <- densidad %>% mutate(densidad = densidad) #factor is the value that you have to use to return to the original values
  return(densidad$densidad)

}
# Maxent_Density <-function(lambda, alpha, minVal=0, maxVal=2, interVals=30, globalZ=1e-3)
# {
#   resultados=NULL
#   
#   x <-seq(minVal,maxVal,length=interVals) #x o S.
# 
#   den <- SME_density(lambda,alpha,x,S=NULL,globalZ)$densidadS
# 
#   return(den)
# }

Maxent_Distribution <- function(lambda, alpha, minVal = 0, maxVal = 2, interVals = 30, globalZ = 1e-3){


  x <- seq(minVal,maxVal,length=interVals) #x o S.
      
  Z         <- SME_density(lambda,alpha, x, S=NULL, globalZ)$Z
  distF     <- SME_distribution(lambda,alpha,x, Z)
  
  
  return(distF)
  
}
  



expectated_value <- function(mat, densi){
  
  x     <- densi[, 1]
  densi <- densi[, -1]
  
  E_x <- NULL
  for (j in 1:nrow(mat)){
    E_x[j] <- trapz(x, x*densi[,j]) 
  }
  
  
  return(mean(E_x, na.rm = T))
}


standard_dev <- function(mat, densi){
  
  x     <- densi[, 1]
  densi <- densi[, -1]
  
  Var_x <- NULL
  for (j in 1:nrow(mat)){
    Var_x[j] <- trapz(x, (x^2)*densi[,j])  - (trapz(x, x*densi[,j]))^2
  }
  
  sd_val <- sqrt(Var_x)
  
  return(mean(sd_val, na.rm = T))
}


