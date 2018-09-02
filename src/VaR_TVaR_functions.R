

##########################################################################
# EMPIRICS

##########################################################################
##############################################
#Empiric VaR   (VaR_function)                #
#Calculate the values of VaR at              #
#gamma% and its confidence intervals at 95%  #                                    
#confidence intervals with                   #
#resampling  
#sample should be greater than 600           #
##############################################

VaR_function<-function(S, gamma, num_rep, sup1=1, inf1=0.80, seed = 1)
{
  set.seed(seed)
  #Calculo del VaR para un valor gamma
  Total=length(S)
  S_order=sort(S, decreasing = FALSE)
  
  #VaR para cada uno de los valores de gamma
  VaR=S_order[round(Total*gamma)]
  
  
  ################################################
  #intervalos de confianza -   remuestreo        #
  ################################################
  
  inf=trunc(inf1*Total)
  sup=trunc(sup1*Total)
  
  #con muestreo sin reemplazo (para sacar intervalos de confianza
  VaR_e=NULL
  for (i in 1:num_rep)
  {
    size=trunc(runif(1,inf,sup)) #muestra entre 10% y 15%
    S_sample=NULL #INICIALIZAR
    S_sample<-sample(S, size, replace = FALSE, prob = NULL)
    S_sample_order<-sort(S_sample, decreasing = FALSE)
    Total<-length(S_sample)
    pos<-Total*(gamma)
    VaR_e[i]<-S_sample_order[round(pos)]
  }
  
  #interval at 95% of confidence
  #x11()
  hist(VaR_e, xlab = 'VaR', ylab = 'frecuencia',main='Histograma de VaR')
  N<-length(VaR_e)
  VaR_e<-sort(VaR_e, decreasing=FALSE)
  #percentiles del 97.5%
  pos_97<-round(N*(0.975))
  #percentiles del 2.5%
  pos_2<-round(N*(0.025))
  abline(v=VaR_e[pos_97], col = "blue")
  abline(v=VaR_e[pos_2], col = "blue", lty=4)
  Int=c(VaR_e[pos_2],VaR_e[pos_97])
  meanVaR=mean(VaR_e)
  
  return(list(VaR=VaR, Interval=Int, meanVaR=meanVaR))
  
} 
#############################################################################

VaR_function_table<-function(S,gamma,sup1=1,inf1=0.80)
{
  if (length(gamma)<1)
  {
    cat(c("gamma should has length grater than zero"),"\n") 
    return(0)
  }
  
  ver=NULL
  ma=matrix(0, nrow = length(gamma), ncol = 7, byrow = FALSE)
  for (i in 1:length(gamma))
  {
    Var_Emp=VaR_function(S[S>0],gamma[i],sup1,inf1)
    TVaR_Emp=TVaR_function(S[S>0],gamma[i],sup1,inf1)
    ver[i]=Var_Emp$VaR # 
    ma[i,1]=gamma[i]
    ma[i,2]=Var_Emp$VaR
    ma[i,3:4]=Var_Emp$Int
    ma[i,5]=TVaR_Emp$TVaR
    ma[i,6:7]=TVaR_Emp$TVaR_Int
  }
  colnames(ma) <- c("gamma","VaR","VaR_inf","VaR_sup","TVaR","TVaR_inf","TVaR_sup")
  ma=data.frame(ma)
  
  return(list(VaR=ver, VaR_table=ma))
}
##############################################################################

VaR_con_function<-function(S,gamma)
{
  if (length(gamma)<1)
  {
    cat(c("gamma should has length greater than zero"),"\n") 
    return(0)
  }
  
  #Calculo del VaR para un valor gamma
  SS=S
  Total=length(SS)
  S_order=sort(SS, decreasing = TRUE)
  
  #VaR para cada uno de los valores de gamma
  VaR=NULL
  for (i in 1:length(gamma))
  {
    VaR[i]=S_order[trunc(Total*(1-gamma[i]))]
  }
  return(VaR)
}
###########################################################################
empirical_distribution_function<-function(S)
{
  Total=length(S)
  S_order=sort(S, decreasing = FALSE)
  F=NULL
  for (i in 1:Total)
  {
    F[i]=length(S_order[S_order<=S_order[i]])/Total
  }
  return(F)
}
############################################################################

VaR_con_function_dist<-function(S,gamma)
{
  if (length(gamma)<1)
  {
    cat(c("gamma should has length grater than zero"),"\n") 
    return(0)
  } 
  
  #Calculo del VaR para un valor gamma
  Total=length(S)
  S_order=sort(S, decreasing = FALSE)
  
  F=NULL
  F=empirical_distribution_function(S)
  FF=NULL
  FF=round(F, digits = 3)
  
  cat("\n")
  cat("\t",c("WARNING: more than one value for VaR"),"\n")
  cat("\n","\t","gamma:","\t","VaR values","\n") 
  #VaR para cada uno de los valores de gamma
  VaR=NULL
  for (i in 1:length(gamma))
  {
    auxi=which(FF==gamma[i])
    if (length(auxi)!=0) 
    {
      V=NULL
      V=S_order[auxi]
      VaR[i]=min(V)
      cat("\t",gamma[i],":","\t",V,"\n") 
    }
    if (length(auxi)==0) 
      VaR[i]=0
  }
  cat("\n")
  return(VaR)
  
}
###########################################################################

#MAXENT VaR & TVaR

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

VaR_function_for_maxent <-function(x,gamma,lambda,alpha,globalZ=0.001){
  # Superior quote of the integral
  maxi <- max(x)*20
  gamma=gamma
  global1=globalZ #0.035 1e-3 #this parameter affects the results (a lot)
  global2=1
  
  y<- seq(from=0,to=1,by=global1)
  f_y<-sum_exp(lambda,alpha,y) #sum of exponentials
  Z<-sintegral(y,f_y)$value #Z normalization parameter
  #Z<-trapz(y,f_y) #Z normalization parameter
  
  x=x[x>0]
  x=sort(x)
  
  U=NULL
  for (i in 1:length(x))
  {
    a=x[i]
    t=seq(a,maxi,global2) 
    f_t=sum_exp(lambda,alpha,exp(-t))
    densidad=exp(-t)*(f_t)*(1/Z)
    integral=sintegral(t,t*densidad-a*densidad)$value
    #integral=trapz(t,t*densidad-a*densidad)
    val=(1/(1-gamma))*(integral)
    U[i]=a+val
  }
  plot(spline(x,U),main = paste("U vs. a (Lemma)"),ylab="U",xlab="x")
  return(list(VaR=x[which.min(U)],TVaR=min(U))) #VaR and TVaR
}
###########################################################################


VaR_function2_poiss<-function(x,gamma,lambda,alpha,maxi=1000,
                              global2=1,N,ita)
{
  
  x=x[x>0]
  x=sort(x)
  
  U=NULL
  for (i in 1:length(x))
  {
    a=x[i]
    t=seq(a,maxi,global2) 
    
    densidad=densityMEMpoisson(lambda,alpha,t,S=NULL,N,ita,meth=5)
    integral=sintegral(t,t*densidad-a*densidad)$value
    #integral=trapz(t,t*densidad-a*densidad)
    val=(1/(1-gamma))*(integral)
    U[i]=a+val
  }
  
  plot(spline(x,U),main = paste("U vs. a (Lemma)"),ylab="U",xlab="x")
  return(list(VaR=x[which.min(U)],TVaR=min(U))) #VaR and TVaR
}


###########################################################################

VaR_function2_exp<-function(x,gamma,lambda,alpha,maxi=1000,
                            global2=1,N,zeta)
{
  
  x=x[x>0]
  x=sort(x)
  
  U=NULL
  for (i in 1:length(x))
  {
    a=x[i]
    t=seq(a,maxi,global2) 
    
    densidad=densityMEM(lambda,alpha,t,S=NULL,N,zeta,mu=NULL,meth=1)$densityMEM
    integral=sintegral(t,t*densidad-a*densidad)$value
    #integral=trapz(t,t*densidad-a*densidad)
    val=(1/(1-gamma))*(integral)
    U[i]=a+val
  }
  
  plot(spline(x,U),main = paste("U vs. a (Lemma)"),ylab="U",xlab="x")
  return(list(VaR=x[which.min(U)],TVaR=min(U))) #VaR and TVaR
}


###########################################################################



#description: find a value F(V)=alpha, using comparations

VaR_maxent<-function(x,gamma,lambda,alpha)
{
  y<- seq(from=0,to=1,by=1e-6) # 3e-3
  sum_e<-source("sum_exp.R")$value #sum of exponentials
  f_y=sum_e(lambda,alpha,y)
  #Z<-sintegral(y,f_y)$value #Z normalization parameter
  Z<-trapz(y,f_y) #Z normalization parameter
  
  
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

nf.method<-function(x,lambda,alpha,gamma,tol,globalZ)
{
  
  fS=SME_density(lambda,alpha,x,S=NULL,globalZ)
  df=fS$densidadS
  Z=fS$Z
  F=SME_distribution(lambda,alpha,x,Z)
  
  func=(F-gamma)^2
  dfunc=2*(F-gamma)*df
  #func=(F-gamma)
  #dfunc=df
  
  if (abs(dfunc)<tol) #10*.Machine$double.eps
  {
    return (x)
  }
  else
  {
    return(x-func/dfunc)
  }
}

###########################################################################

#Return the value F(V)=gamma using Newton_Raphson method
#where 
# g(x)= (F(V)-gamma)^2
# g'(x))= 2(F(V)-gamma)F'(V)
# or
# g(x)= (F(V)-gamma)
# g'(x))= F'(V)

Newton_Raphson_min<-function(lambda,alpha,gamma,initial_value,tol=1e-5,iter=1000,globalZ=0.001)
{
  tam=length(initial_value)
  if (tam!=1)
  {
    cat(c("We need only a value, no more"),"\n") 
    return(0)
  }
  
  start<-initial_value
  col=rainbow(20)
  xn<-NULL
  xn<-c(xn,start) 
  
  n=1
  
  while (iter>0)
  {
    n=n+1
    xn<-c(xn,nf.method(xn[n-1],lambda,alpha,gamma,tol,globalZ))
    if (xn[n]<0)
    {
      cat(c("Try a new initial value"),"\n")
      return(0)
    }
    if(abs(xn[n]-xn[n-1])<tol) break #100*.Machine$double.eps
    iter=iter-1
  }
  if (iter>0)
    return(xn[length(xn)])
  if (iter==0)
    return(0)
  
}

############################################################################

#############################################################################

#TVar functions

##############################################
#Empiric TVaR   (TVaR_function)               #
#Calculate the values of VaR at              #
#gamma% and its confidence intervals at 95%  #
#option 1: confidence intervals with the     #
#density                                     #
#option 2: confidence intervals with         #
#resampling                                  #
##############################################

TVaR_function<-function(S,gamma,sup1=1,inf1=0.80)
{
  
  #TVaR calculation
  SS=S
  Total=length(SS)
  S_order=sort(SS, decreasing = FALSE)
  
  p2=trunc(Total*gamma)
  p1=1/(Total-p2+1)
  
  suma=0
  for (j in (p2):Total)
    suma=suma+S_order[j]
  
  
  TVAR=p1*suma
  
  ######################
  #confidence interval
  inf=trunc(inf1*Total)
  sup=trunc(sup1*Total)
  TVAR_values=NULL
  for (i in 1:1000)
  {
    
    size=trunc(runif(1,inf,sup)) #muestra entre 90% y 100%
    sev=NULL #INICIALIZAR
    sev<-sample(SS, size, replace = TRUE, prob = NULL)
    
    #ordenar los datos generados
    Sev<-sort(sev, decreasing=FALSE)
    Total<-length(Sev)
    pos<-Total*(gamma)
    
    #para el TVar
    
    p2=round(Total*(gamma)) #trunc(Total*(gamma))
    p1=1/(Total-p2+1)
    
    suma=0
    for (j in (p2):Total)
      suma=suma+Sev[j]
    
    TVAR_values[i]=p1*suma
    
  }
  
  #histograma TVaR
  #x11()
  hist(TVAR_values,xlab = 'TVaR', ylab = 'frecuencia',main='Histogram de TVaR')
  N<-length(TVAR_values)
  TVAR_values<-sort(TVAR_values, decreasing=FALSE)
  #percentiles del 97.5%
  pos_97<-round(N*(0.975))
  #percentiles del 2.5%
  pos_2<-round(N*(0.025))
  abline(v=TVAR_values[pos_97], col = "blue")
  abline(v=TVAR_values[pos_2], col = "blue", lty=4)
  Int=c(TVAR_values[pos_2],TVAR_values[pos_97])
  
  return(list(TVaR=TVAR,TVaR_Int=Int))
  
}