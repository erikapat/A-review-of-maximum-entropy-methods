
#objetive function
#return a scalar

function(lambda,alpha,mu,trap=1)
{

 global<<-0.0001
 #verifications
 # lambda,alpha size bigger than 0
 if (length(alpha)<=0 | length(lambda)<=0 | length(mu)<=0) return(1001)
 if (length(alpha)!=length(lambda) | length(alpha)!=length(mu) | length(mu)!=length(lambda)) return(1001)

 #empiezan los calculos
 #the objective function have two parts: an integral and escalar product
 #integrate density from 0 to 1
 estimate=NULL
 
 #y<-seq(0,1,length=100)
 y<- seq(from=0,to=1,by=global) #*

 #call the sum_of exponential
 #f<-source("sum_exp.R")$value
 fy<-sum_exp(lambda,alpha,y)
 #fy=(fy[fy!=Inf]) #para evitar que se vuelva muy grande
 if (length(fy)==1) if (fy==1001) return(1001)

 if (trap==0)
 estimate<-sintegral(y,fy)$value #Z 
 if (trap==1)
 estimate<-trapz(y,fy) 
 
 #objetive function
 #original
 #norma<-source("norma.R")$value
 estimate2<-log(estimate)+sum(lambda*mu)#+((0.001)/2)*(norma(lambda))^2 
 #print(sum(lambda*mu))
 return(estimate2)
}

#verifications
#the size of lambda, mu and alpha should be equal.