#objetive function FOR THE METHOD OF MAXIMUN ENTROPY IN THE MEAN
#FUNCTION:funMEM

function(lambda, alpha, mu, zeta = 10, N = 100){

 Z  <- Zmem(lambda, alpha, mu, N, zeta)
 if (Z == -999){
   return(-999)
}
 if (Z < 0){
  cat("Change constant zeta, log(Z) = NaN, Z ", Z, "\n")
  return(-999)
 }

 estimate <-  log(Z) + sum(lambda*mu) 

 return(estimate)
}

#verifications
#the size of lambda, mu and alpha should be equal.