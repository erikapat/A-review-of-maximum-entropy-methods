
# extra values of SMEE

# #################
# # density
# #################
# #density and graphs of S and estimated
# source("SME_functions.r")
# x=seq(0, max(S) + 0.1 ,length= 100) #x o S.
# den=SME_density(lambda,alpha,x,S,globalZ=1e-4)
# Z=den$Z
# fS=den$densidadS
# 
# sum_e<-source("sum_exp.R")$value #sum of exponentials
# y<- seq(from=0,to=1,0.0015)
# f_y=sum_e(lambda,alpha,y)
# mus=NULL
# Z<-trapz(y,f_y)
# for (i in 1:length(lambda))
# {
#   mus[i]<-trapz(y,(y^alpha[i])*(f_y)/Z)
# }
# mus
# mp$int_mu
# 
# pk=exp(-mu_interval[,1]*lambda)/(exp(-mu_interval[,1]*lambda)+exp(-mu_interval[,2]*lambda))
# ek=(mu_interval[,1]*exp(-mu_interval[,1]*lambda)+mu_interval[,2]*exp(-mu_interval[,2]*lambda))/(exp(-mu_interval[,1]*lambda)+exp(-mu_interval[,2]*lambda))
# pk
# ek
# 
# cuenta=mus+pk*mu_interval[,1]+(1-pk)*mu_interval[,2]
# cuenta

#---------------------------------------------------------------------------------------------------------------------------------------------------------
theorical_interval_mu<-function(S1,k=8,M=20,rep=1000,alph=0.05)
{
  
  #calculation of alpha
  dec=1.5
  dec1=rep(dec,k)
  kk1=1:k
  alpha=dec1/kk1
  
  mu_matrix=matrix(0,ncol=k,nrow=rep)
  for (j in 1:rep)
  {
    set.seed(j)
    indS=NULL
    indS=sample(1:length(S1),size=M)
    S=NULL
    S=S1[indS]
    
    #Prob of N=0
    po=length(S[S==0])/length(S)
    laplace=rep(0,k)
    mu=rep(0,k)
    #laplace transform
    for (i in 1:k)
      laplace[i]=(1/length(S))*sum(exp(-alpha[i]*S)) 
    #moments
    mu=(laplace-po)/(1-po)
    
    mu_matrix[j,]=mu
  }
  
  means=NULL
  sds=NULL
  for (i in 1:k)
  {
    means[i]=mean(mu_matrix[,i])
    sds[i]=sd(mu_matrix[,i])
  }
  
  #mu_inf=means-qnorm(1-(alph)/2)*(sds/sqrt(rep))
  #mu_sup=means+qnorm(1-(alph)/2)*(sds/sqrt(rep))
  mu_inf=-qnorm(1-(alph)/2)*(sds/sqrt(rep))
  mu_sup=qnorm(1-(alph)/2)*(sds/sqrt(rep))
  int_mu=cbind(mu_inf,mu_sup)
  
  return(int_mu)
}



#others (check)
###########################################################################
theorical_interval_mu_xx<-function(S1,k=8,M=20,rep=1000,alph=0.05)
{
  #calculation of alpha
  dec=1.5
  dec1=rep(dec,k)
  kk1=1:k
  alpha=dec1/kk1
  
  mu_matrix=matrix(0,ncol=k,nrow=rep)
  for (j in 1:rep)
  {
    set.seed(j)
    indS=NULL
    indS=sample(1:length(S1),size=M)
    S=NULL
    S=S1[indS]
    
    #Prob of N=0
    po=length(S[S==0])/length(S)
    laplace=rep(0,k)
    mu=rep(0,k)
    #laplace transform
    for (i in 1:k)
      laplace[i]=(1/length(S))*sum(exp(-alpha[i]*S)) 
    #moments
    mu=(laplace-po)/(1-po)
    mu_matrix[j,]=mu
  }
  
  means=NULL
  sds=NULL
  for (i in 1:k)
  {
    means[i]=mean(mu_matrix[,i])
    sds[i]=sd(mu_matrix[,i])
  }
  
  #mu_inf=means-qnorm(1-(alph)/2)*(sds/sqrt(rep))
  #mu_sup=means+qnorm(1-(alph)/2)*(sds/sqrt(rep))
  mu_inf=-qnorm(1-(alph)/2)*(sds/sqrt(rep))
  mu_sup=qnorm(1-(alph)/2)*(sds/sqrt(rep))
  int_mu=cbind(mu_inf,mu_sup)
  
  return(int_mu)
}


#--------------------------------------------------------------------------------
# Modificado enero 2016 (resto del codigo del 2015)
theoretical_interval_mu_new <-function(S1,k=8,M=20,rep=1000,alph=0.05){
  #calculation of alpha
  dec=1.5
  dec1=rep(dec,k)
  kk1=1:k
  alpha=dec1/kk1
  
  mu_matrix <- matrix(0,ncol = k, nrow = rep)
  for (j in 1:rep){
    set.seed(j)
    indS <- NULL
    indS <- sample(1:length(S1), size = M)
    S <- NULL
    S <- S1[indS]
    
    #Prob of N = 0
    po      <- length( S[S==0] )/length(S)
    laplace <- rep(0, k)
    mu      <- rep(0, k)
    #laplace transform
    for (i in 1:k)
      laplace[i] <- (1/length(S))*sum(exp(-alpha[i]*S)) 
    #moments
    mu <- (laplace-po)/(1-po)
    mu_matrix[j,] <- mu
  }
  
  means <- NULL
  sds   <- NULL
  for (i in 1:k)
  {
    means[i] <- mean(mu_matrix[,i])
    sds[i]   <- sd(mu_matrix[,i])
  }
  
  #mu_inf=means-qnorm(1-(alph)/2)*(sds/sqrt(rep))
  #mu_sup=means+qnorm(1-(alph)/2)*(sds/sqrt(rep))
  mu_inf <- -qnorm(1-(alph)/2)*(sds/sqrt(rep))
  mu_sup <-  qnorm(1-(alph)/2)*(sds/sqrt(rep))
  #mu_inf <-   abs(mu_data - means) + qnorm((1 - alph)/2)*(sds/sqrt(rep))
  #mu_sup <-   abs(mu_data - means) - qnorm((1 - alph)/2)*(sds/sqrt(rep))
  
  int_mu <-  cbind(mu_inf, mu_sup)
  
  return(int_mu)
}