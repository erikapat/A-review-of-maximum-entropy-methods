



#general functions:
#source("general_functions.r")

library(pracma)

#################################################################
#DESCRIPTION:                                                   #
#Simulation of the compound distribution of losses              #
#from a Poisson and a log-normal                                #
#Parameters:                                                    #
#INPUT                                                          #
#name: name of the file where you save the data                 #
#ele: parameter of the poisson                                  #
#mu and sigma: lognormal parameters                             #
#M: number of values generated from the Poisson                 #
#OUTPUT:                                                        #
#file with the generated data                                   #
#################################################################

simul_compound<-function(name, ele, me, sigma, M, seed = 1)
{
  #simulation of Poiss-lognormal
  #M data from a poisson of lambda =1, k=1,...M
  M_ori = M
  M = 2*M

  set.seed(seed)
  #dist Poisson
  n=rpois(M,ele) #N
  
  S=vector(mode = "numeric", length = M)
  indi=NULL
  # for each k 
  for (k in 1:M){
    xc = (rlnorm(n[k], meanlog = me, sdlog = sigma)) 
    S[k]=sum(xc)
    indi=c(indi,xc)
    #meanlog, sdlog: mean and standard deviation of the distribution on the log scale with default values of 0 and 1 respectively.
  }
  
  n <- n[S>0]
  S <- S[S>0]
  indi <- indi[indi>0]
  index <- sample(1:length(S), M_ori)
  S    <- S[index]
  indi <- indi[index]
  n    <- n[index]
  #GRAPH of Individuals
  hist(indi,breaks=15,freq=FALSE,include.lowest = TRUE, main =  'Individuals - logNormal distribution')
  
  #GRAPH of the simulation
  hist(S, breaks=15, freq=FALSE, include.lowest = TRUE, main =  'Compounded distributions')
  h <- hist(S, breaks=100, plot = F)
  true.densi <- data.frame(x = h$mids, densidad = h$density)
  #save
  write.table(S, file = paste0('data/', name, 'compounded_dist.dat'), row.names = FALSE, col.names = FALSE)
  write.table(indi, file = paste0("data/", name, ".individal_lognormal.dat"), row.names = FALSE, col.names = FALSE)
  write.table(true.densi, file = paste0("data/", name, "true.densi.dat"), row.names = FALSE, col.names = FALSE)
  write.table(n, file = paste0("data/", name, "frequency.dat"), row.names = FALSE, col.names = FALSE)
  
  return(0) 
}
#########################################################################


simul_compound_gamma<-function(name, ele, a, b, M, seed = 1)
{
  #simulation of Poiss-lognormal
  #M data from a poisson of lambda =1, k=1,...M
  
  #dist Poisson
  set.seed(seed)
  n=rpois(M,ele) #N
  S=vector(mode = "numeric", length = M)
  indi=NULL
  # for each k 
  for (k in 1:M)
  {
    xc=(rgamma(n[k], shape = a, rate = b)) 
    S[k]=sum(xc)
    indi=c(indi,xc)
  }
  
  # #GRAPH of Individuals
  # x11()
  # hist(indi[indi>0],breaks=15,freq=FALSE,include.lowest = TRUE)
  # 
  # #GRAPH of the simulation
  # x11()
  # hist(S[S>0],breaks=15,freq=FALSE,include.lowest = TRUE)
  
  #save
  # write.table(S, file = name,row.names = FALSE,col.names = FALSE)
  # write.table(indi, file = "indi2.dat",row.names = FALSE,col.names = FALSE)
  #save
  write.table(S, file = paste0('data/', name, 'compounded_dist.dat'), row.names = FALSE, col.names = FALSE)
  write.table(indi, file = paste0("data/", name, ".individual_GAMMA.dat"), row.names = FALSE, col.names = FALSE)
  #write.table(true.densi, file = paste0("data/", name, "true.densi.dat"), row.names = FALSE, col.names = FALSE)
  write.table(n, file = paste0("data/", name, "frequency.dat"), row.names = FALSE, col.names = FALSE)
  return(0) 
}

#########################################################################
########################################################################

#the same simulation like before, with the diference that the parameters 
#of the distributions are always changing.

simul_compound_dif<-function(name, M)
{
  #simulation of Poiss-lognormal
  #M data from a poisson of lambda =1, k=1,...M
  
  #dist Poisson
  ele=runif(M,1,5)
  n=rpois(M,ele) #N
  me=runif(M,0,0.5)
  sigma=runif(M,0.5,0.75)
  
  S=vector(mode = "numeric", length = M)
  
  # for each k 
  for (k in 1:M)
  {
    xc=(rlnorm(n[k], meanlog = me[k], sdlog = sigma[k])) 
    S[k]=sum(xc)
  }
  
  #GRAPH of the simulation
  hist(S[S>0],breaks=10,freq=FALSE,include.lowest = TRUE)
  
  #save
  write.table(S, file = name,row.names = FALSE,col.names = FALSE)
  
  return(0) 
}

##########################################################################

#HESSIAN MATRIX

hessian_matrix<-function(alpha,lambda,x,mu_interval)
{
  #HESSIAN FOR THE METHODS (2.1) OF THE PAPER: Density reconstruction with errors in the data.
  
  #Calculation of Z
  sum_e<-source("sum_exp.R")$value #sum of exponentials
  y<- seq(from=0,to=1,by=1e-4)
  f_y=sum_e(lambda,alpha,y)
  #Z<-sintegral(y,f_y)$value #Z normalization parameter
  Z<-trapz(y,f_y)
  
  #Calculation of the final density
  #Calculation of the final density
  f_y=sum_e(lambda,alpha,y)
  densidad=(f_y)*(1/Z)
  
  H=matrix(0,ncol=length(alpha),nrow=length(alpha))
  for (i in 1:length(alpha))
  {
    for (j in 1:length(alpha))
    {
      term1=(y^(alpha[i]+alpha[j]))*densidad
      term2=(y^(alpha[i]))*densidad
      term3=(y^(alpha[j]))*densidad
      H[i,j]=trapz(y,term1)-trapz(y,term2)*trapz(y,term3)
    }
  }
  
  return(H)
  
}

###############################################################################################
################################################################################################

hessian_matrix_2.2<-function(alpha,lambda,pk,mu_interval)
{
  #HESSIAN FOR THE METHODS (2.2) OF THE PAPER: Density reconstruction with errors in the data.
  
  #Calculation of Z
  sum_e<-source("sum_exp.R")$value #sum of exponentials
  y<- seq(from=0,to=1,by=1e-4)
  f_y=sum_e(lambda,alpha,y)
  #Z<-sintegral(y,f_y)$value #Z normalization parameter
  Z<-trapz(y,f_y)
  
  #Calculation of the final density
  f_y=sum_e(lambda,alpha,y)
  densidad=(f_y)*(1/Z)
  
  H=matrix(0,ncol=length(alpha),nrow=length(alpha))
  for (i in 1:length(alpha))
  {
    for (j in 1:length(alpha))
    {
      term1=(y^(alpha[i]+alpha[j]))*densidad
      term2=(y^(alpha[i]))*densidad
      term3=(y^(alpha[j]))*densidad
      H[i,j]=trapz(y,term1)-trapz(y,term2)*trapz(y,term3)
      if(i==j)
        H[i,j]=H[i,j]+(mu_interval[i,1]^2)*pk[i]*(1-pk[i])+(mu_interval[i,2]^2)*pk[i]*(1-pk[i])-2*(mu_interval[i,2]*mu_interval[i,1])*pk[i]*(1-pk[i])
    }
  }
  
  return(H)
  
}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#EJEMPLO 1 
#simulacion de una pois-lognormal  ---> paper xxxx: log-normal: meanlog = 0.5, sdlog = 0.75, set.seed(5000), N = 5000
                                                                #pois  : lambda 1

por_lognormal_simulations <- function(N = 5000, mean_log = .5, sd_log = .75){
  options(digits= 16)

  set.seed(5000)
  
  S=NULL
  k=1
  while ( k<5000 )
  {
    n <- rpois(1,1)
    if (n!=0)
    {
      xc   <- (rlnorm(n, meanlog = mean_log, sdlog = sd_log))
      S[k] <- sum(xc)
      k    <- k+1
    }
  }
  
  hist(S)
  min(S)
  max(S)
  hist(S, breaks = 10,freq  =FALSE, include.lowest = TRUE)
}
