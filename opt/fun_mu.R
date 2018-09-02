#objetive function
#return a scalar

function(lambda,alpha,mu,trap=1)
{

 global=0.0001

 #empiezan los calculos

 #it is necessary install the function
 #install.packages("Bolstad", lib="/my/own/R-packages/")
 library(Bolstad)
 library(pracma)
 #the objective function have two parts: an integral and escalar product
 #integrate density from 0 to 1
 estimate=NULL
 
 y<- seq(from=0,to=1,by=global) #*

 #call the sum_of exponential
 #f<-source("sum_exp.R")$value
 fy<-sum_exp(lambda,alpha,y)
 if (length(fy)==1) if (fy==1001) return(1001)

 if (trap==0)
   estimate<-sintegral(y,fy)$value #Z 
 if (trap==1)
   estimate<-trapz(y,fy) 
 
 #objetive function
 #norma<-source("norma.R")$value

 est1=sum((mu[,2]-mu[,1])*abs(lambda)/2)
 est2=-sum((mu[,2]+mu[,1])*lambda/2)
 estimate2<-log(estimate)+est1-est2

 #print(estimate2)
 return(estimate2)
}

#verifications
#the size of lambda, mu and alpha should be equal.