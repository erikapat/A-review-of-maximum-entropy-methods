#POISSON & EXPONENTIAL FUNCTIONS -------------------------------------------------------------------------

# POISSON FUNCTION

# FUNCTIONS TO USE BEFORE DE OPTIMIZATION PROCESS
#----------------------------------------------------------------------------------------------------

#MATRIX A: 
matrixA <- function(alpha,N)
{
  k=length(alpha)
  A=matrix(0, nrow=k, ncol=N)
  for (i in 1:k) 
  {
    for (j in 1:N)
      A[i,j]=((2*j-1)/(2*N))^alpha[i] #(1/N)*
  }
  return(A)
}


#------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# FUNCTIONS TO USE AFTER DE OPTIMIZATION PROCESS

densityMEMpoisson <- function(lambda, alpha, x, N = 200, ita = 2, meth = 5){
  options(warn = -1) #dont show Warnings 0 or a negative number
  if (meth < 1 | meth > 6) 
  {
    cat( "Error: Invalid interpolation method" , "\n")
    return(0)
  }
  
  A     <- matrixA(alpha,N)
  #Maxentropic estimator of the solution 
  #(it is not yet the solution until the change of variables) 
  z  <- sum(exp(-t(A)%*%lambda))
  if (z == 0){
    cat("\n", "\t","Error: NaN results", "\n")
    return(0)
  }
  xOpt <- exp(-t(A)%*%lambda)/z   #z!=Z
  if (sum(xOpt)!=1) cat("\n","WARNING: sum of x* not equal to one ","\n") #the sum should be equal to 1
  #print(sum(xOpt))
  xOpt <- N*xOpt
  
  
  #INTERPOLACION
  j    <- seq(1, length(xOpt))
  yAux <- (j - 1)/N
  
  #Cubic Hermite spline
  metho <- c("fmm", "periodic", "natural", "linear", "constant")
  
  xOpt_aux <- sort.int(xOpt, index.return = TRUE)
  
  
  if (meth == 6){
    i <- aspline(yAux[xOpt_aux$ix], xOpt_aux$x, xout = exp(-x), method = "improved", degree = 100) #method: "original" or "improved"
  } else if (meth == 4 | meth == 5){
    i <- approx(yAux[xOpt_aux$ix], xOpt_aux$x, xout = exp(-x), method = metho[meth]) #method: "linear" or "constant"
  } else if (meth==1 | meth== 2 | meth==3){
    i <- spline(yAux[xOpt_aux$ix], xOpt_aux$x,xout = exp(-x), method = metho[meth], ties = mean) #"periodic" gives a better approximation but gives Warnings
  }                                                  #other options that does not gives so good results is "fmm"
  densidad <- exp(-x)*i$y
  
  densidad <- data.frame(x = x, densidad = densidad)
  return(densidadPoisson = densidad)
}

#FITTED DISTRIBUTION FUNCTION
MEMpoisson_distribution<-function(x,lambda,alpha,N,ita,meth=1)
{
  x=sort(x)
  
  #Calculation of the distribution function
  integral=NULL
  for (i in 1:length(x))
  {
    t=seq(0,x[i],0.001) 
    densidad= densityMEMpoisson(lambda, alpha, t, N = 200, ita = 2, meth = 5)$densidad
    integral[i]=sintegral(t, densidad)$value #distribution function
    #WARNING: if I use trapz the result does not exist.
    #integral[i]=trapz(t,densidad) #distribution function
  }
  F_maxent = data.frame(x = x, F = integral)
  return(F_maxent)
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# EXPONENTIAL FUNCTION

#FUNCTION: Zmem
Zmem <- function(lambda,alpha,mu,N,zeta){
  #we star here
  #values for some values:
  k=length(alpha)
  
  A     <- matrixA(alpha,N)
  
  Z=1
  for (j in 1:N)
  {
    aux=zeta/(zeta+A[,j]%*%lambda)
    Z=Z*aux
  }
  
  return(Z)
}

#----



###################################################################
#densityMEM_exp:                                                   #
#Define the values of the exponential-MEM                          #
#                                                                  #
#INPUT                                                             #
# lambda: optimum lambda                                           #
# alpha:  alpha vector                                             #
# x    : values of the x-axis  
#N 
#zeta=10, 
#mu = NULL, 
#meth = 1#

#OUTPUT                                                            #
# Z : normalization factor                                         #
# densidadS: density                                               #
####################################################################

densityMEM_exp <- function(lambda, alpha, x, N = 1, zeta = 1, mu = NULL, meth = 4)
{
  k  <- length(lambda)
  Zm <- source("Zmem.R")$value
  Z  <- Zm(lambda, alpha, mu, N, zeta)
  mA <- source("matrixA.R")$value
  A  <- mA(alpha, N)
  
  #Maxentropic estimator of the solution 
  #(it is not yet the solution until the change of variables) 
  xOpt  <- (1/(zeta + (t(A)%*%lambda))) #(eq(12))
  #print(sum(xOpt))
  
  error <- "NA"
  if (length(mu) > 0){
    error <- norma(A%*%xOpt - mu) 
  }
  #INTERPOLACION
  j    <- seq(1, length(xOpt))
  yAux <- (j - 1)/N
  
  #Cubic Hermite spline
  metho <- c("fmm", "periodic", "natural", "linear", "constant")
  
  library(akima)
  
  if (meth == 6){
    i <- aspline(yAux, N*xOpt, xout = exp( -x ), method = "improved", degree = 100) #method: "original" or "improved"
  } else if (meth == 4 | meth == 5){
    i <- approx(yAux, N*xOpt, xout = exp(-x), method = metho[meth]) #method: "linear" or "constant"
  } else if (meth == 1 | meth== 2 | meth==3){
    i <- spline(yAux, N*xOpt, xout = exp(-x), method = metho[meth], ties = mean) #"periodic" gives a better approximation but gives Warnings
  }                                                   #other options that does not gives so good results is "fmm"
  
  densidad <- exp(-x)*i$y
  
  return(list(densityMEM = densidad,errorMEM = error))
}


