
#gradient of the objetive function
#a vector with k elements with the value of the gradient
#Name: gradientMEM2 (POISSON)

function(lambda, alpha, mu, correction, ita=2, N = 200) {
 #norma <- source("norma.R")$value
 #verifications
 # lambda,alpha size bigger than 0
 if (length(alpha) <=0 | length(lambda) <= 0 | length(mu) <= 0) {
   return(1000)
}else if (length(alpha) != length(lambda) | length(alpha) != length(mu) | length(mu) != length(lambda)) {
  return(1001)
}
 #We start here
 #...for some values...
 k  <- length(alpha)
 A     <- matrixA(alpha,N)
 
 xOpt <- exp(-t(A)%*%lambda)
 aux  <- A%*%xOpt
 aux2 <- -ita*aux
 grad <- rep(0,k)

 if (correction == 0)
 {
   for (i in 1:k){
     grad[i]= aux2[i,1] + mu[i] 
   }
 } else {
   for (i in 1:k){
     grad[i] <- aux2[i,1] + mu[i]  + (1e-3)*lambda[i] 
   }
 }

 return(grad)
}

#verifications
#the size of lambda, mu and alpha should be equal.