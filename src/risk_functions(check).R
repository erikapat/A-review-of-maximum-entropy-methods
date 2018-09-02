
###########################################################################
#description: find a value F(V)=alpha, using compararisons

VaR_maxent<-function(x,gamma,lambda,alpha)
{
  y<- seq(from=0,to=1,by=1e-6) # 3e-3
  sum_e<-source("sum_exp.R")$value #sum of exponentials
  f_y=sum_e(lambda,alpha,y)
  Z<-sintegral(y,f_y)$value #Z normalization parameter
  #Z<-trapz(y,f_y) #Z normalization parameter
  
  
  #Calculation of the distribution function
  x=x[x>0]
  x=sort(x,decreasing = TRUE)
  int=NULL
  integral=NULL
  for (i in 1:length(x))
  {
    t=seq(0,x[i],0.001) 
    f_t=sum_e(lambda,alpha,exp(-t)) 
    densidad=exp(-t)*(f_t)*(1/Z)
    integral=trapz(t,densidad) #distribution function
    int[i]=round(integral, digits = 4)
  }
  
  gamma=signif(gamma,3)
  final_VaR=NULL
  auxi=NULL
  for (j in 1:length(gamma))
  {
    auxi=which(int==gamma[j])
    if (length(auxi)!=0)
      final_VaR[j]=min(x[which(int==gamma[j])]) #vector of indexes
    if (length(auxi)==0)
      final_VaR[j]=0 #not answer
    print(x[which(int==gamma[j])])
  }
  return(final_VaR)
} #NOTE: for values of gamma>=.99 this function is not good enough

###########################################################################

#MAXENTROPICOS

##########################################################################

#VaR_function2
#####################################################################
#Calculation of VaR (second version)                                #
#This use the lemma 5.1 in the paper                                #
#Gzyl H., Tagliani A. Determination of the distribution of total    #
#loss from the fractional moments of its exponential.               #
#Applied Mathematics and Computation 219 (2012): 2124-2133.         #
#####################################################################
#maxi=cota superior de la integral

VaR_function_theorem<-function(x,gamma,lambda,alpha,maxi=1000,globalZ=0.001,global2=0.1)
{
  global1=globalZ #0.035 1e-3 #this parameter affects the results (a lot)
  
  y<- seq(from=0,to=1,by=global1)
  sum_e<-source("sum_exp.R")$value #sum of exponentials
  f_y=sum_e(lambda,alpha,y)
  #Z<-sintegral(y,f_y)$value #Z normalization parameter
  Z<-trapz(y,f_y) #Z normalization parameter
  x=x[x>0]
  x=sort(x)
  
  U=rep(0,length(x))
  for (i in 1:length(x))
  {
    a        <- x[i]
    t        <- seq(a,maxi,global2) 
    f_t      <- sum_e(lambda,alpha,exp(-t))
    densidad <- exp(-t)*(f_t)*(1/Z)
    densidad <- densidad[!is.na(densidad)]
    t        <- t[!is.na(densidad)]
      #integral=sintegral(t,t*densidad-a*densidad)$value
      integral=trapz(t,t*densidad-a*densidad)
      val=(1/(1-gamma))*(integral)
      U[i]=a+val

  }
  
  plot(x,U,main = paste("U vs. a (Lemma)"),ylab="U",xlab="x")
  return(list(VaR=x[which.min(U)],TVaR=min(U))) #VaR and TVaR
}

#---------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------

VaR_function_theorem2<-function(x,gamma,lambda,alpha,maxi=1000,globalZ=0.001,global2=0.1)
{
  
  
  global1=globalZ #0.035 1e-3 #this parameter affects the results (a lot)
  
  y<- seq(from=0,to=1,by=global1)
  sum_e<-source("sum_exp.R")$value #sum of exponentials
  f_y=sum_e(lambda,alpha,y)
  #Z<-sintegral(y,f_y)$value #Z normalization parameter
  Z<-trapz(y,f_y) #Z normalization parameter
  x=x[x>0]
  x=sort(x)
  
  U <- mapply(U_curve, 1:length(x))

  plot(x,U,main = paste("U vs. a (Lemma)"),ylab="U",xlab="x")
  return(list(VaR=x[which.min(U)],TVaR=min(U))) #VaR and TVaR
}

#-------------------------------------------------------------------------------------------------------------------------------------

U_curve <- function(i){
  a=x[i]
  t=seq(a,maxi,global2) 
  f_t=sum_e(lambda,alpha,exp(-t))
  densidad=exp(-t)*(f_t)*(1/Z)
  densidad <- densidad[!is.na(densidad)]
  t        <- t[!is.na(densidad)]
  #integral=sintegral(t,t*densidad-a*densidad)$value
  integral=trapz(t,t*densidad-a*densidad)
  val=(1/(1-gamma))*(integral)
  U=a+val
  
  return(U)
}



