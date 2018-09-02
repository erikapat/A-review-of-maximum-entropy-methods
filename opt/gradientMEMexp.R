#gradient of the objetive function
#a vector with k elements with the value of the gradient
#Name: gradientMEM

function(lambda, alpha, mu, correction, zeta = 10, N = 100) {
 
 #norma <- source("norma.R")$value
 k     <- length(alpha)

 A     <- matrixA(alpha,N)
 
 xOpt  <- 1/(zeta + t(A)%*%lambda)
 aux   <- A%*%xOpt
 grad  <- rep(0,k)
 for (i in 1:k){
   grad[i] <- -aux[i,1] + mu[i]  
 }
 if (correction == 1){
   grad <- grad + (1e-3)*lambda
 }
 return(grad)
}

#verifications
#the size of lambda, mu and alpha should be equal.