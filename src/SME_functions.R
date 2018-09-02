


#FUnctions that are required for the optimization procedure (before) ---------------------------------------------------------

sum_exp <- function(lambda,alpha,y) {
  
  #sum of exponencials, needed for the gradient and objetive function
  #lambda,alpha: vector of the same dimension
  #if all is correct the return is a vector, if not the return is a scalar (this means a error)
  
  #verifications
  # lambda,alpha size bigger than 0
  if (length(alpha)<=0 | length(lambda)<=0 | length(y)<=0) return(1001)
  
  if (length(alpha)!=length(lambda)) return(1001)
  
  
  #empiezan los calculos
  #y<-seq(0,1,length=100)
  sumfy<-rep(0,length(y))
  for (j in 1:length(alpha)) 
  {
    sfy<-lambda[j]*(y^alpha[j])
    sumfy<-sumfy+sfy
    sfy<-NULL
  }
  
  f<-exp(-sumfy) 
  
  return(f)
  
}


# Functions that use the result of the optimization method to generates densities, CDF's and validation's figures.

SME_density<-function(lambda, alpha, x, globalZ=1e-3){
  # Calculate the density using optimal lambda
  # lambda and alpha are the parameters used to get the density
  # points use to calculate its density
  # globalZ: this improves the results in the most simple cases, and the the more difficult cases help the results without error
  #(usually you play with this parameter)
   y   <- seq(from=0,to=1, by = globalZ)
   f_y <- sum_exp(lambda,alpha,y)
  
  #Z normalization parameter
  Z<-trapz(y,f_y)
  if (is.na(Z)==TRUE){
    cat("Change variable GlobalZ","\n")
    stop()
    geterrmessage() 
    return(0)
  }
  
  #Calculation of the final density
  f_x      <- sum_exp(lambda,alpha,exp(-x))
  densidad <- exp(-x)*(f_x)*(1/Z)
  
  densidadS <- data.frame(x = x, densidad = densidad)
  
  #OUTPUT
  return(list(Z = Z, densidadS=densidadS))
  
}

#------------------------------------------------------------------------------------------------------------------------

#FITTED DISTRIBUTION FUNCTION
#Z :       Normalization parameter
# lambda: value obtained from the optimization problem
# alpha: defined before the optimization problem
# x:     points of the densities to be obtained
SME_distribution <- function(lambda, alpha, x, Z, num_points_x = 100)
{
  cat('[Info] Calculating F(x) for each element of x')
  #Calculation of the distribution function
  integral <- NULL
  x = sort(x)
  #apply map
  integral <- mapply(simple_dist, 1:length(x), MoreArgs = list(lambda, alpha, x, Z, num_points_x))
  integral[integral>=1] <- 1 #greater than 1
  return(integral)
}

simple_dist <- function(i, lambda, alpha, x, Z, num_points_x){
  
  t           <- seq(0, x[i], length = num_points_x) 
  f_t         <- sum_exp(lambda,alpha,exp(-t)) 
  densidad    <- exp(-t)*(f_t)*(1/Z)
  return(trapz(t, densidad)) #distribution function
}

#-------------------------------------------------------------------------------------------------------------------------

# distribution function of F
#S :      datos...
# lambda: value obtained from the optimization problem
# alpha: defined before the optimization problem

F_diff <- function(S, lambda, alpha, globalZ = 0.0055){
  
  options(digits = 4)
  SS <- sort(S)

  # get the density
  jj  <- SME_density(lambda,alpha,SS,globalZ)
  F   <- SME_distribution(lambda, alpha, SS, Z=jj$Z)
  
  F   <- round(F, digits = 3) #fitted distribution

  F_maxent = data.frame(x = SS, F = F)
  
  return(list(F_maxent = F_maxent))
  
}

#--------------------------------------------------------------------------
empirical_CDF <- function(S){
  options(digits = 4)
  SS <- sort(S)
  
  kl <- seq(1,length(SS),length = length(SS)) 
  CDF <- (length(SS)-kl+0.5)/length(SS)
  CDF <- rev(CDF)
  F_e <- CDF
  F_e <- round(F_e,digits = 3) #empiric distribution
  F_true = data.frame(x = SS, F = F_e)
  
  return(list(F_true = F_true))
}

#------------------------------------------------------------------------
