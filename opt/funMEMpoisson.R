#objetive function FOR THE METHOD OF MAXIMUN ENTROPY IN THE MEAN
#FUNCTION:funMEM2

function(lambda,alpha,mu,ita=2,N=200) {

 #verifications
 # lambda,alpha size bigger than 0
 if (length(alpha) <= 0 | length(lambda) <= 0 | length(mu) <= 0) return(1001)
 if (length(alpha) != length(lambda) | length(alpha) != length(mu) | length(mu) != length(lambda)) return(1001)

 #Zm<-source("Zmem.r")$value
 #Z=Zm(lambda,alpha,mu,N,zeta)
 #if (Z==-999)
 #  return(-999)

 A     <- matrixA(alpha,N)

 suma1 <- 1 - exp(-t(A)%*%lambda)
 suma2 <- sum(suma1)
 suma3 <- -ita*(suma2) + sum(lambda*mu) #+((1e-3)/2)*(norma(lambda))^2

 #Z=prod(exp(-ita*(1-exp(-t(A)%*%lambda)))) 
 #print(Z) 
 #estimate=log(Z)+sum(lambda*mu)
 #print(estimate)

 return(suma3)
}

#verifications
#the size of lambda, mu and alpha should be equal.