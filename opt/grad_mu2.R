
#gradient of the objetive function
#a vector with k elements with the value of the gradient

function(lambda,alpha,mu,correction,trap=1)
{
 #it is necessary install the function
 #install.packages("Bolstad", lib="/my/own/R-packages/")

 global=0.0001 
 mu_int=read.table(file = "rds/mu_matrix.dat")
 k=length(lambda)
 #in this case \mu is a matrix
 if (dim(mu_int)[1]!=k | dim(mu_int)[2]!=2 ) return(1001)

 #verifications
 # lambda,alpha size bigger than 0
 if (length(alpha)<=0 | length(lambda)<=0) return(1000)

 if (length(alpha)!=length(lambda)) return(1001)
 #the objective function have two parts: an integral and escalar product
 #integrate density from 0 to 1
 estimate=NULL
 y<- seq(from=0,to=1,by=global) #*

 #value of Z
 #call the sum_of exponential
 #sum_e<-source("sum_exp.R")$value
 noZ<-sum_exp(lambda,alpha,y)
 if (length(noZ)==1) if (noZ==1001) return(1001)

 if (trap==0)
    Z<-sintegral(y,noZ)$value #Z 
 if (trap==1)
    Z<-trapz(y,noZ) 

 #exp of the sum
 sum<-sum_exp(lambda,alpha,y)
 if (length(sum)==1) if (sum==1001) return(1001)

 estimate<-NULL
 for (i in 1:length(alpha))
 {
  fy<-(-y^alpha[i])*sum
  if (trap==0)
    int_fy=sintegral(y,fy)$value
  if (trap==1)
    int_fy<-trapz(y,fy)

 
  if (correction==0){
      estimate[i]<-(1/Z)*int_fy+mu[i]-(mu_int[i,1]*exp(-mu_int[i,1]*lambda[i])+ mu_int[i,2]*exp(-mu_int[i,2]*lambda[i]))/(exp(-mu_int[i,1]*lambda[i])+ exp(-mu_int[i,2]*lambda[i]))
      #estimate[i]<-(1/Z)*int_fy+mu[i]-(mu_int[i,1]*exp(-mu_int[i,1]*lambda[i])+ mu_int[i,2]*exp(-mu_int[i,2]*lambda[i])) 
    }
  if (correction==1){
     #estimate[i]<-(1/Z)*int_fy+mu[i]-(mu_int[i,1]*exp(-mu_int[i,1]*lambda[i])+ mu_int[i,2]*exp(-mu_int[i,2]*lambda[i])) +0.001*lambda[i]
      estimate[i]<-(1/Z)*int_fy+mu[i]-(mu_int[i,1]*exp(-mu_int[i,1]*lambda[i])+ mu_int[i,2]*exp(-mu_int[i,2]*lambda[i]))/(exp(-mu_int[i,1]*lambda[i])+ exp(-mu_int[i,2]*lambda[i]))
      +0.001*lambda[i]
  }
 }

 #print(norma(estimate))
 return(estimate)
}

#verifications
#the size of lambda, mu and alpha should be equal.
