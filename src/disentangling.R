


#source("disentanglingRr")
#INPUT
#xx1, xx2: random variables
#hh1, hh2: frequencies of the random variables

#(xx2,hh2): random variables and frerquencies together
#(xx1,hh1):

library(actuar)

dis_data2<-function(xx2, hh2, xx1, hh1) 
{ 
 x<-matrix(0, nrow =max(c(length(xx1),length(xx2))),ncol=2, byrow=FALSE)


 #EMPIEZO A TRABAJAR CON EL VECTOR MAS PEQUENO
 x1=NULL
 X2=NULL
 cc=NULL
 hh=NULL
 if (length(xx2)<=length(xx1))
 {
   x2=xx2
   hh=hh2
   x1=xx1
   cc=hh1
 }
 if (length(xx2)>=length(xx1))
 {
   x2=xx1
   hh=hh1
   x1=xx2
   cc=hh2
 }

x[1:length(x2),1]=x2 #el vector de indices mas pequenios primero.
x[1:length(x2),2]=hh



ult=NULL
ult2=NULL
for (i in 1:length(x2))
{
 for (j in 1:length(x1))
 {
  if (x1[j]==x2[i])
  {
   x[i,1]=x2[i]
   x[i,2]=cc[j]+hh[i]
   ult=i
   ult2=j
  }
 }
}

if (is.null(ult))
{
 ult=length(x2)
 ult2=0
}


sec=seq((ult2+1),length(x1))


x[(ult+1):(ult+length(sec)),1]=x1[(ult2+1):length(x1)]
x[(ult+1):(ult+length(sec)),2]=cc[(ult2+1):length(x1)]

index=which(x[,1]!=0)
xfinal=x[1:max(index),]

return(xfinal)


}

#########################################################################################

dis_data<-function(n1,n2) 
{ 
 l=0
 kn=c(n1,n2)
 mat=matrix(0, nrow =(max(kn)+1),ncol=2, byrow=FALSE)
 for (i in 1:(max(kn)+1))
 {
   mat[i,2]=length(kn[kn==l])
   mat[i,1]=l
   l=l+1
 }
 return(mat)
}
##############################################################################
#########################################################################################

dis_data_improve<-function(kn)   #kn=c(n1,n2) 
{ 
 l=0
 #kn=c(n1,n2)
 mat=matrix(0, nrow =(max(kn)+1),ncol=2, byrow=FALSE)
 for (i in 1:(max(kn)+1))
 {
   mat[i,2]=length(kn[kn==l])
   mat[i,1]=l
   l=l+1
 }
 return(mat)
}



#########################################################################################

recursion<-function(x,limy=NULL)
{

 col=x[,2]

 coln=col[col!=0]
 indice=x[,1]
 indicen=indice[col!=0]

 akb=rep(0,length(coln)-1)

 for (i in 1:length(coln)-1)
 {
  akb[i]=(coln[i+1]/coln[i])*indicen[i+1]
 }

 kk=akb
 #akb=round(akb)
 akb=c(0,t(t(akb)))

 #par(mfrow=c(1,2))
 if (length(limy)==0)
    plot(indicen[-1],akb[2:length(akb)],xlab = 'k', ylab = 'kr(k) ')
 if (length(limy)!=0)
    plot(indicen[-1],akb[2:length(akb)],xlab = 'k', ylab = 'kr(k) ',ylim=c(0,limy))

 matr1 <- matrix(c(indicen,coln,akb), nrow = length(akb), ncol=3, byrow=FALSE)

 return(matr1)

}

########################################################################################
bin_parameters<-function(matr1)
{
   y1=matr1[,2]
   x1=matr1[,1]
   recta1=lm(y1~x1)

   coef1=recta1$coefficients

   a_coef=coef1[1]
   b_coef=coef1[2]

  beta=b_coef/(1-b_coef)
  r=(a_coef*(1+beta))/beta+1
  p0=(1+beta)^(-r)
  media=r*beta
  vari=r*beta*(1+beta)

  x11()
  plot(x1,y1)
  lines(x1,a_coef+b_coef*x1)
  return(list(beta=beta, r=r,p0=p0,media=media, vari=vari))

}

#######################################################################################
#######################
#      METODO EM      #
#######################



#matr1: recursion matrix
#mod: model --> "EII", "VII", "EEI", "EVI", "VEI", "VVI"
#G: number of groups



EM_method<-function(matr1,mod=2,G=2)
{

 if (mod<=0) 
 {
    cat(c("Error: Invalid Model Name "),"\n")
    return(0)
 }
 if (mod>6) 
 {
    cat(c("Error: Invalid Model Name "),"\n")
    return(0)
 }
    
 models<-c("EII", "VII", "EEI", "EVI", "VEI", "VVI")

 library(mclust)
 a=Mclust(matr1[,-2],G, modelNames=models[mod]) #"VII"
 x11()
 plot(a,matr1[,-2],what="BIC")

 #calculation of probabilities...
 cl1=NULL
 Negentropy=NULL
 ncol=dim(matr1)[1]
 if (G==1)
 {
  Negentropy=1
  Entropy=0
 }
 if (G>1)
 {
  cl1=a$z
 
 # print(cl1)
 #cat("\n","Prob:","\n") #A matrix whose [i,k]th entry is the probability that observation i in the test data belongs to the kth class.
 #print(round(cl1,2))
 
  suma=0
  for (i in 1:ncol)
  {
   for (j in 1:G)
   {
    suma=suma+(-cl1[i,j]*log(cl1[i,j]))
   }
  }
  Negentropy=1-suma/(ncol*log(G))
  Entropy=suma
 }
 print(Negentropy)
 

 a$bic
 a$BIC
 a$loglik
 a$parameters
 table(a$class)

 return(list(lista=a$class, BIC=a$BIC,bic=a$bic,
 loglikelihood=a$loglik,probas=cl1,degrees_of_freedom=a$df, entropy=Entropy,
  negentropy=Negentropy))
}

##########################################################################################

#CRITERIOS PARA DETERMINAR EL N?MERO DE GRUPOS EN ALGORITMOS DE AGLOMERAMIENTO
#..............................................
#VALORES DE ENTRADA:
#log_lik: log verosimilitud
#grupos: n?mero de grupos a considerar
#n: total de datos
#ENT: entropia
#printi= 1 si quieres imprimir en pantalla, 1 es el valor por defecto
#..............................................
#VALORES DE SALIDA
#result= matriz con los resultados de los criterios AIC, BIC y ICL-BIC

num_grupos<-function(log_lik,ENT, grupos,n,degrees_of_freedom,printi=1) 
{
 #d n?mero de coef a estimar
    #d=4*grupos-1
    d=degrees_of_freedom
    AIC=-2*log_lik+2*d
    AIC3=-2*log_lik+3*d
    BIC=-2*log_lik+d*log(n)
    #ENT %entropia
    ICL_BIC=BIC+2*ENT
 #end

names=c('  AIC  ' , '  AIC3 ', '  BIC  ', 'ICL_BIC');
result=c(AIC, AIC3, BIC, ICL_BIC);

if (printi==1)
{
    cat(names,'\n')
    cat(result,'\n')
}

}
########################################################################################


#matr1: recursion matrix
#mod: model --> "EII", "VII", "EEI", "EVI", "VEI", "VVI"
#G: number of groups
#dist: type of distribution
#lines: drawn the regression lines
#kmeans= 0 no k-means analysis
#        1    k-means analysis
#ylim2: max value of y axis

do_groups<-function(matr1,list,dist="Negative Binomial",lineas=1,kmeans=0,ylim2=NULL,true0=0,xlabeli='k',ylabeli='kr(k)',ylim1=0){
  beta=NULL
  r=NULL
  p0=NULL
  media=NULL
  vari=NULL
  n=NULL
  p=NULL

 if (length(ylim2)==0){
  ylim2=max(matr1[-1,3])*2
 }
 xini=1
 if (true0==1){
  xini=0
 }

 #plot(matr1[-1,1],matr1[-1,3])
 plot(matr1[-1,1],matr1[-1,3], type = "n",xlab = xlabeli, ylab = ylabeli,ylim=c(ylim1,ylim2),xlim=c(xini,max(matr1[-1,1])))  # setting up coord.  system #,xlim=c(1,max(matr1[,1])+2)

 x=matr1[,1]
 y=matr1[,3]

 a_coef=NULL #intercepto a
 b_coef=NULL # b
 
 G=max(list)
 for (i in 1:G){
  points(matr1[list==i,1],matr1[list==i,3],pch=5*i)
  if (length(matr1[list==i,1])>1)
  {
  
   x1=x[list==i]
   y1=y[list==i]

   recta1=lm(y1~x1)

   coef1=recta1$coefficients

   a_coef[i]=coef1[1]
   b_coef[i]=coef1[2]

    if (lineas==1)
     lines(x1,a_coef[i]+b_coef[i]*x1)
  }

  if (length(matr1[list==i,1])==1){
   a_coef[i]=matr1[list==i,3] #intercepto
   b_coef[i]=0
  }
 }

 
 if (dist=="Negative Binomial"){
  beta=b_coef/(1-b_coef)
  r=(a_coef*(1+beta))/beta+1
  p0=(1+beta)^(-r)
  media=r*beta
  vari=media*(1+beta)  
 }

 if (dist=="Poisson"){
  media=a_coef
  vari=a_coef
  p0=exp(-media)
 }
 if (dist=="Binomial"){
  p=-b_coef/(1-b_coef)
  n=(a_coef/(-b_coef))-1
  media=n*p
  vari=n*p*(1-p)
  p0=(1-p)^n
 }

 if (kmeans==1){
 #######################
 #      K-MEANS        #
 #######################

 #data=matrix(c(x,y), nrow=length(x), ncol=2, byrow=F)

 clusters=kmeans(matr1[,-2],G,iter.max=100,nstart=25)
 points(clusters$centers, pch = c(8,4), cex=5)

 cat(c("############################################################### "),"\n")
 cat("\n")
 cat(c("K-means: "),"\n")
 cat(c("Centers: ", clusters$centers,"\n"))
 cat("\n")
 cat(c("############################################################### "),"\n")

 print(aggregate(matr1[,-2],by=list(clusters$cluster),FUN=mean))
}



 return(list(a_coef=a_coef, b_coef=b_coef, beta=beta, r=r,p0=p0,n=n, p=p,media=media, vari=vari))
}

##########################################################################################
#calculation of proyeccions with a rotation matrix (I don't use this)

proyecciones<-function(x, alpha)
{
 dime=dim(x)
 
 xx=NULL
 yy=NULL
 for (i in 1:dime[1])
 {
   xx[i]=(x[i,1]*cos(alpha)+x[i,2]*sin(alpha))
   yy[i]=(-x[i,1]*sin(alpha)+x[i,2]*cos(alpha))
 }

 h=rep(0,dime[1])
 par(mfrow = c(1, 2))
 plot(xx,yy,pch=5) 
 plot(h,yy,pch=10)

 return(0)
}

################################################################################

#L1 and L2 norm between the real and estimated density
#L1 and L2 norm divided by n

density_errors1<-function(density_real, density_est)
{
   
  diff=sum(abs(density_real-density_est))/length(density_real)
  diff2=sqrt(sum((density_real-density_est)^2))/length(density_real)
  
  return(list(L1=diff, L2=diff2))

}

##################################################################################


density_simple_errors_SME<-function(lambda,alpha,S,globalZ=1e-3,bins=15,ejemplo=5)
{
 #S es respecto a quien comparas (data empirica)
 #globalZ:sensibilidad de Z
 #Calculation of Z

 sep=bins #number of bins

 d=hist(S[S>0],breaks=sep,include.lowest = TRUE,
           plot=FALSE) #main = paste("20 bars")

 densidad_simulada=d$density

 #Z
 sum_e<-source("sum_exp.R")$value #sum of exponentials
 y<- seq(from=0,to=1,by=globalZ)
 f_y=sum_e(lambda,alpha,y)
 #Z<-sintegral(y,f_y)$value #Z normalization parameter
 Z<-trapz(y,f_y)


 int=d$breaks
 print(int)
 int[1]=int[1]+0.01 #OJO PARA EVITAR NAN's en el caso de la pareto u otras dist.
 
 #error real vs. maxentropic 
 ll1=NULL
 densidad1=NULL
 densidad_real1=NULL
 aux1=0
 aux11=0
 for ( i in 1:(length(int)-1) )
 {
  ll=seq(int[i],int[i+1],0.01)

  f_ll=sum_e(lambda,alpha,exp(-ll))
  densidad=exp(-ll)*(f_ll)*(1/Z)
  set.seed(1)
  ###########################################################
  #densidad_real=dlnorm(ll, meanlog=par1, sdlog=par2) #****CAMBIAR
  if (ejemplo==5)
  {
   densidad_real=(10/16)*dpareto(ll, shape=5, scale=1, log = FALSE)+(6/16)*dpareto(ll,    shape=15,   scale=2.5, log = FALSE)
  }
  if (ejemplo==3)
  {
   densidad_real=2/10*dlnorm(ll, meanlog = -1, sdlog = 0.25, log = FALSE) +8/10*dlnorm(ll,    meanlog = -1, sdlog =0.25, log = FALSE)
  }
    if (ejemplo==4)
  {
  densidad_real=(4/19)*dpareto(ll, shape=15, scale=2, log = FALSE)+(15/19)*dpareto(ll,   shape=15, scale=2, log = FALSE)
  }
    if (ejemplo==6)
  {
   densidad_real=4/10*dgamma(ll,shape=5,rate=15) +6/10*dgamma(ll,shape=5,rate=15)
  }
    if (ejemplo==1)
  {
   densidad_real=1/3*dfrechet(ll,  scale =2, shape=4)+2/3*dbeta(ll, shape1=1, shape2=5) 
  }
    if (ejemplo==2)
  {
   mu2=7.5
   mu3=19.5
   mu4=27
   mu5=32
   mu=mu2+mu3+mu4+mu5
   densidad_real=(mu2/mu)*dbeta(ll, shape1=1, shape2=25) + (mu3/mu)*dweibull(ll, shape=1, scale = 0.1)  +(mu4/mu)*dfrechet(ll, scale =0.01, shape=2)+(mu5/mu)*dgamma(ll, shape=0.1, rate = 3)

 }
  ###########################################################

  #aux1=aux1+trapz(ll,abs(densidad_real-densidad))
  diff=abs(densidad_real-densidad)
  aux1=aux1+sintegral(ll,diff)$value
  diff2=(densidad_real-densidad)^2
  aux11=aux11+sintegral(ll,diff2)$value
  #aux11=aux11+trapz(ll,(densidad_real-densidad)^2)
  

  ll1=c(ll1,ll)
  densidad1=c(densidad1,densidad)
  densidad_real1=c(densidad_real1,densidad_real)
 }
  
    x11()
   d=hist(S[S>0],breaks=sep,freq=FALSE,include.lowest = TRUE,
          ylim = c(0,max(c(densidad,d$density,densidad_real))+0.1), xlim=c(0, max(c(S[S>0],x))),
          xlab="S", main = paste("")) #main = paste("20 bars")
   lines(spline(ll1,densidad1),col="black",lwd =3,lty=2) 
   lines(spline(ll1,densidad_real1),col="red",lwd =3,lty=2) 
   legend(max(c(S[S>0],x[densidad>0]))-5,max(c(densidad,d$density))-0.02, c("SME","REAL"), lty=c(2,2), lwd=c(2,2),col=c("black","red")) 


  

 
 #error simulada vs. maxentropic 
 aux2=0
 aux22=0
 aux3=0
 aux33=0
 for ( i in 1:(length(int)-1) )
 {
  if (d$counts[i]<2) next
  if (i<(length(int))) sup=int[i+1]
  
  val=S[S<sup]
  val=val[val>int[i]]
  densi= rep(densidad_simulada[i],length(val))


  f_val=sum_e(lambda,alpha,exp(-val))
  densidad=exp(-val)*(f_val)*(1/Z)
  set.seed(1)
  ##############################################################
  #densidad_real=dlnorm(val, meanlog=par1, sdlog=par2) #******CAMBIAR
  if (ejemplo==5)
  {
  densidad_real=10/16*dpareto(val, shape=5, scale=1, log = FALSE)+6/16*dpareto(val,   shape=15, scale=2.5, log = FALSE)
  }
  if (ejemplo==3)
  {
     densidad_real=2/10*dlnorm(val, meanlog = -1, sdlog = 0.25, log = FALSE) +8/10*dlnorm(val,  meanlog = -1, sdlog =0.25, log = FALSE)
  }
    if (ejemplo==4)
  {
  densidad_real=(4/19)*dpareto(val, shape=15, scale=2, log = FALSE)+(15/19)*dpareto(val,   shape=15, scale=2, log = FALSE)
  }
    if (ejemplo==6)
  {
   densidad_real=(4/10)*dgamma(val,shape=5,rate=15) +(6/10)*dgamma(val,shape=5,rate=15)
  }
    if (ejemplo==1)
  {
   densidad_real=1/3*dfrechet(val,  scale =2, shape=4)+2/3*dbeta(val, shape1=1, shape2=5) 
  }
    if (ejemplo==2)
  {
   mu2=7.5
   mu3=19.5
   mu4=27
   mu5=32
   mu=mu2+mu3+mu4+mu5
   densidad_real=(mu2/mu)*dbeta(val, shape1=1, shape2=25) + (mu3/mu)*dweibull(val, shape=1, scale = 0.1)  +(mu4/mu)*dfrechet(val, scale =0.01, shape=2)+(mu5/mu)*dgamma(val, shape=0.1, rate = 3)

 }
  ##############################################################

  #sintegral(val,densidad)$value
  #aux2=aux2+trapz(val,abs(densi-densidad))
  dif=abs(densi-densidad)
  aux2=aux2+sintegral(val,dif)$value
  dif2=(densi-densidad)^2
  aux22=aux22+sintegral(val,dif2)$value

  difh=abs(densi-densidad_real)

  
  aux3=aux3+sintegral(val,difh)$value
  difh2=(densi-densidad_real)^2
  aux33=aux33+sintegral(val,difh2)$value
  

 }

 # ll=seq(int[length(int)],int[length(int)]*10,0.01)
 #Calculation of the final density
 #f_x=sum_e(lambda,alpha,exp(-ll))
 #densidad=exp(-ll)*(f_x)*(1/Z)
 # aux2=aux2+trapz(ll,abs(densidad))
 # aux22=aux22+trapz(ll,(densidad)^2)

   error2=aux2
   error22=sqrt(aux22)
 
   error3=aux3
   error33=sqrt(aux33)

  error1=aux1 
  error11=sqrt(aux11)


   #OUTPUT
   return(list(Z=Z,L1_real=error1,L1_sim=error2,
               L2_real=error11,L2_sim=error22,  L1_real_sim=error3,
               L2_real_sim=error33))

  #L1_real & L2_real: maxent density & real density 
  #L1_sim  & L2_sim : maxent density & simulated density 
  #L1_real_sim & L2_real_sim: simulated density & real density 

}

###########################################################################################



density_simple_errors_SME2<-function(S,lambda,alpha,gamma,globalZ=0.0055,ejemplo=5)
{

 
 SS=sort(S[S>0.01])

 kl=seq(1,length(SS),length=length(SS)) 
 CDF=(length(SS)-kl+0.5)/length(SS)
 CDF=rev(CDF)
 F_e=CDF

 ##########################################################################
 jj=density_simple_errors_SME(lambda,alpha,S=SS,globalZ=1e-3,bins=15,ejemplo)
 ##########################################################################
 
 F=SME_distribution_l(lambda,alpha,SS,Z=jj$Z)
 set.seed(1)
 ##############################################################################
 #F_real=plnorm(SS, meanlog = me, sdlog =sd, lower.tail = TRUE, log.p = FALSE)
 if (ejemplo==5)
 F_real=(10/16)*ppareto(SS, shape=5, scale=1, log = FALSE)+(6/16)*ppareto(SS, shape=15, scale=2.5, log = FALSE)
 if (ejemplo==3)
   F_real=2/10*plnorm(SS, meanlog = -1, sdlog = 0.25, log = FALSE) +8/10*plnorm(SS,    meanlog = -1, sdlog =0.25, log = FALSE)
    if (ejemplo==4)
  {
  F_real=(4/19)*ppareto(SS, shape=15, scale=2, log = FALSE)+(15/19)*ppareto(SS,   shape=15, scale=2, log = FALSE)
  }
    if (ejemplo==6)
  {
   F_real=4/10*pgamma(SS,shape=5,rate=15) +6/10*pgamma(SS,shape=5,rate=15)
  }
    if (ejemplo==1)
  {
   F_real=1/3*pfrechet(SS,  scale =2, shape=4)+2/3*pbeta(SS, shape1=1, shape2=5) 
  }
    if (ejemplo==2)
  {
   mu2=7.67
   mu3=19.2
   mu4=27.3
   mu5=32.37
   mu=mu2+mu3+mu4+mu5
   F_real=(mu2/mu)*pbeta(SS, shape1=1, shape2=25) + (mu3/mu)*pweibull(SS, shape=1, scale = 0.1)  +(mu4/mu)*pfrechet(SS, scale =0.01, shape=2)+(mu5/mu)*pgamma(SS, shape=0.1, rate = 3)

 }

 ##############################################################################
 x11()
 plot(SS,F_e,type="l", lwd = 3.5,lty=1,main="Empirical and SME distribution", ylab="F",xlab="S")
 lines(SS,F,col="gray",lty=2, lwd = 3.5)
 #lines(SS,F_real,col="red",lty=2, lwd = 3.5)
 legend(max(SS)-1.2,0.2, c("Empirical","SME"), lty=c(1,4), lwd=c(3.5,2.5),col=c("black","gray")) 


 #MEAN ABSOLUTE ERROR #L1-norm
 MAE_sim=(1/length(F))*sum(abs(F-F_e))
 
 #L2-norm
 #ROOT MEAN SQUARED ERROR (#BRIER SCORE: measures the mean squared error. Negative orientation (smaller score the better))
 MSE_sim=sqrt((1/length(F))*sum((F-F_e)^2))

 #MEAN ABSOLUTE ERROR #L1-norm
 MAE_real=(1/length(F))*sum(abs(F-F_real))
 
 #L2-norm
 #ROOT MEAN SQUARED ERROR (#BRIER SCORE: measures the mean squared error. Negative orientation (smaller score the better))
 MSE_real=sqrt((1/length(F))*sum((F-F_real)^2))

 #MEAN ABSOLUTE ERROR #L1-norm
 MAE_real_sim=(1/length(F))*sum(abs(F_e-F_real))
 
 #L2-norm
 #ROOT MEAN SQUARED ERROR (#BRIER SCORE: measures the mean squared error. Negative orientation (smaller score the better))
 MSE_real_sim=sqrt((1/length(F))*sum((F_e-F_real)^2))


 return(list(F_emp=F_e,F_fit=F,MAE_sim=MAE_sim,RMSE_sim=MSE_sim,MAE_real=MAE_real,
             RMSE_real=MSE_real,MAE_real_sim=MAE_real_sim,
             RMSE_real_sim=MSE_real_sim))

}
####################################################################################


SME_distribution_l<-function(lambda,alpha,x,Z)
{
 sum_e<-source("sum_exp.R")$value #sum of exponentials

 #Calculation of the distribution function
 integral=NULL
 for (i in 1:length(x))
 {
  t=seq(0,x[i],0.001) 
  f_t=sum_e(lambda,alpha,exp(-t)) 
  densidad=exp(-t)*(f_t)*(1/Z)
  #integral[i]=sintegral(t,densidad)$value #distribution function
  integral[i]=trapz(t,densidad) #distribution function
 }
 integral[integral>=1]=1 #greater than 1
 return(integral)
}

#################################################################################
#SIMULATION OF MIXED SEVERITIES

#sizeS=size of the sample
#a and b, interval of the simulation

mix_simul<-function(sizeS=500, a=0, b=2,rep=10)
{
 library(actuar)
 sim=NULL
 for (j in 1:sizeS)
 {
  x=NULL
  s=1
  while (s<rep)
  {
  
   r1=runif(1,0,1)
   r2=runif(1,0,1)
   c=a+(b-a)*r1

  #puede cambiarse, introduciendo un archivo
  ########################################################################################
   den=(10/16)*ppareto(c, shape=5, scale=1, log = FALSE)+(6/16)*ppareto(c, shape=15, scale=2.5, log = FALSE)
#den=0.99*pgamma(c,100,1/80)+0.01*pgamma(c,100*100,1/80)
  ########################################################################################

  if (den>=r2)
  {
   x[s]=c
   s=s+1
  }
   
  } 
 
 
  sim[j]=min(x)
 }
  hist(sim,freq = FALSE,xlab='X', ylab='Frequency', main=" Individual Losess")
  return(sim)
}

#########################################################################################
#SIMULATION OF MIXED SEVERITIES

#sizeS=size of the sample
#a and b, interval of the simulation

mix_simul2<-function(sizeS=500, seed=1)
{
 library(actuar)
 set.seed(seed)
 sim=NULL

   uu <- runif(sizeS)   
   den1=rpareto(sizeS, shape=5, scale=1)
   den2=rpareto(sizeS, shape=15,scale=2.5)
   sim<- ifelse(uu < 10/16, den1,den2)
  hist(sim,freq = FALSE,xlab='X', ylab='Frequency', main=" Individual Losess")
  return(sim)
}

##########################################################################################

#SIMULATION OF MIXED SEVERITIES

#sizeS=size of the sample
#a and b, interval of the simulation

mix_simul2_log<-function(sizeS=500, seed=1)
{
 library(actuar)
 set.seed(seed)
 sim=NULL

   uu <- runif(sizeS)   
   den1=rlnorm(sizeS, meanlog = -1, sdlog = 0.25)
   den2=rlnorm(sizeS, meanlog = -1, sdlog = 0.25)
   sim<- ifelse(uu < 2/8, den1,den2)
  hist(sim,freq = FALSE,xlab='X', ylab='Frequency', main=" Individual Losess")
  return(sim)
}


#########################################################################################
#sizeS=size of the sample
#a and b, interval of the simulation

mix_simul_gamma<-function(size1=500,size2=400, seed=1)
{
 library(actuar)
 set.seed(seed)
 sim=NULL
 sizeS=min(size1,size2)
 
   uu <- runif(sizeS)   
   den1=rgamma(sizeS,shape=5,rate=15)
   den2=rgamma(sizeS,shape=5,rate=15)
   sim<- ifelse(uu < 4/19, den1,den2)
   cen=0
   adi=length(sim)
   while (cen!=abs(size1-size2))
   {
       adi=adi+1
       cen=cen+1
       uu <- runif(1) 
       den2=rgamma(1,shape=5,rate=15)  
       sim[adi]=ifelse(uu < 4/19, 0,den2)
       if (sim[adi]==0) 
       {
        adi=adi-1
        cent=cen-1
       }
   }

   print(length(sim))
  hist(sim,freq = FALSE,xlab='X', ylab='Frequency', main=" Individual Losess")
  return(sim)
}

#########################################################################################
#SIMULATION OF MIXED SEVERITIES

#sizeS=size of the sample
#a and b, interval of the simulation

mix_simul_pareto_2<-function(sizeS=500, seed=1)
{
 library(actuar)
 set.seed(seed)
 sim=NULL

   uu <- runif(sizeS)   
   den1=rpareto(sizeS, shape=15, scale=2)
   den2=rpareto(sizeS, shape=15,scale=2)
   sim<- ifelse(uu < 4/19, den1,den2)
  hist(sim,freq = FALSE,xlab='X', ylab='Frequency', main=" Individual Losess")
  return(sim)
}

##############################################################################
#Boostrap for the calculation of mu, calculation of interval for mu (laplace transform de individual losses)
#source("SME_functions.R")

#gamma=0.95 #interval of 95% of confidence

boot_mu<-function(S,k,iter=1000,subsample=100,replacement=0,gamma=0.95,globalZ=1e-2)
{

 M=length(S)
 #calculation of alpha
 dec=1.5
 dec1=rep(dec,k)
 kk1=1:k
 alpha=dec1/kk1

 #calculation of moments mu
 
 

  mu=rep(0,k)
  hu=NULL
 for (i in 1:iter)
 {
  laplace=rep(0,k)

  sub_S=sample(S, size=subsample, replace = replacement, prob = NULL)
  #laplace transform
  for (i in 1:k)
  {
   laplace[i]=(1/length(sub_S))*sum(exp(-alpha[i]*sub_S)) 
   #Prob of N=0
   po=length(sub_S[sub_S==0])/length(sub_S)
  }

  aux=(laplace-po)/(1-po)
  
  

  ############
  #OJO CAMBIAR
  ############
 # lval=laplace_individuals(lambda=NULL,alpha=NULL,distribution="Poisson",Lambda=10,globalZ=globalZ,mu=aux)

 lval=laplace_individuals(lambda=NULL,alpha=NULL,distribution="Binomial",prob=0.76,size=114,globalZ=1e-7,mu=aux)

  #lval=laplace_individuals(lambda=NULL,alpha=NULL,distribution="negative binomial",prob=0.4270,size=11.77,globalZ=1e-2,mu=aux)
  hu=cbind(hu,lval)
  #moments
  mu=mu+aux
  
 }

 did=((1-gamma)/2)
 nn=dim(hu)[2]

 
 ifi=trunc(did*nn)

 sp=trunc((1-did)*nn)
 fil=NULL
 inter=NULL
 hu_m=NULL
 for (i in 1:dim(hu)[1])
 {
   fil=hu[i,]
   fil=sort(fil)
   hu_m[i]=mean(fil)
   vc=NULL
   vc=c(fil[ifi],fil[sp])
   inter=rbind(inter,vc)
 }


return(list(inter=inter,hu_m=hu_m))

}








