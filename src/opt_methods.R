
#opt_met.R: contains the following functions:
#BB():
#error():
#Rosembrock():
#maxent_parameters:


#################################################################
#DESCRIPTION:                                                   #
#BB-METHOD                                                      #
#Parameters:                                                    #
#INPUT                                                          #
#lambda:initial guess of values for lambda                      #
#alpha                                                          #
#mu: moments                                                    #
#alpha, mu, lambda have the same size                           #
#M: maximum number of iterations (default 1000)                 #
#tolerance: Stopping criterion (default 1e-5)                   #
#OUTPUT:                                                        #
#lambda_Actual=optimum lambda                                   #
# or errors: -999 (NaN values)                                  #
#            -666 (number of iterations reach)                  #
#             000 (infinite loop)                               #
#################################################################
#ultimo cambio los sigma.


BB<-function(lambda,alpha,mu,nameFun="fun.R",nameGrad="gradient.R",M=1000,tolerance=1e-5,step=1,correction=0)
{
 start.time <- Sys.time() #starting to measuring time
 options(digits=8)
 #functions to use:
 #call the functions
 gradient<-source(nameGrad)$value
 f<-source(nameFun)$value
 #norma<-source("norma.r")$value
 #delta<-source("delta.r")$value

 #Initiate variables:
 lambda_Actual=lambda
 grad=gradient(lambda,alpha,mu,correction=correction)
 grado=norma(grad)
 iter=1
 funcEval=f(lambda,alpha,mu)
 gamma=1e-4
 epsilon=1e-10
 sigma1=0.1 #0.1 #***
 sigma2=0.5 #0.5
 beta=1
 N=10 #***
 betas=NULL #***
 m=19 #***
 
 for (i in 1:length(grad))  
 if (is.na(grad[i])==T) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-999))

 dif=1
 dif_grad=1
 #STEP 1:
 while(norma(grad)>tolerance & iter<M)
 #while(dif_grad>tolerance & iter<M) #stop criteria
 {
  #STEP 2:
  if (beta<=epsilon | beta>=1/epsilon) beta=delta(grad)

  #STEP 3:
  phi=1/beta
 
  #STEP 4: Nonmonotone line search to find the right lambda
  #sk=-grad search direction
  #phi= stepsize, depends of the method used 
  lambda_aux=lambda_Actual-phi*grad #iteration formula
  
  fE=f(lambda_aux,alpha,mu)

  cen=0
  while(is.na(fE)==T)
  {
   if (cen==20) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-999))
   cen=cen+1
   #"Error: lambdas out of range"
   lambda_aux=lambda_Actual-(0.5^(cen))*grad
   fE=f(lambda_aux,alpha,mu)
  } #the idea with this loop is control lambda

  if (fE==-999) 
  {
    #print("Error: lambda out of the range")
    return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=000))
    #negative values (funMEM), infinite loop
  } 

  #descending method
  index=iter-(0:min(iter-1,N))
  maximo=max(funcEval[index])-gamma*(grad%*%grad)*(phi) #*(-phi) 
  cen=0
  while(fE>maximo) #descend condition
  {
   #STEP 5:
   cen=cen+1
   set.seed(100)
   sigma=runif(1)*(sigma2-sigma1)+sigma1
   phi=sigma*phi
   lambda_aux=lambda_Actual-phi*grad
   fE=f(lambda_aux,alpha,mu)
   if (is.na(fE)==T) return(list(lambda_Opt=0,iter=iter,grad=norma(grad),error=-999))
   if (cen==100) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=000))
  }

  #STEP 6: update beta
  lambda_old=lambda_Actual
  lambda_Actual=lambda_aux
  yk=gradient(lambda_Actual,alpha,mu,correction=correction)-grad

  #the step is controlled by this...
  if (step>4 | step<1)
  {
    cat("Error: Not valid step","\n")
    return(list(lambda_Opt=0,iter=iter,grad=norma(grad),error=999.9))
  }
  if (step==1)
  beta=-grad%*%yk/(phi*grad%*%grad)
  if (step==2)
  {
   sk=-grad
   beta=yk%*%yk/(phi*sk%*%yk)
  } 
  if (step==3)
  {
   sk=-grad
   beta2=yk%*%yk/(phi*sk%*%yk)
   beta1=-grad%*%yk/(phi*grad%*%grad)
   beta=beta1
   if(beta1/beta2<=0.25) beta=beta2
  } 
  if (step==4)
  {
   sk=-grad
   beta2=yk%*%yk/(phi*sk%*%yk)
   betas[iter]=1/beta2
   sl=seq(max(1,iter-m),iter)
   beta1=min(betas[sl])
   beta=1/beta1
  } 
  if (is.na(beta)==TRUE)
  {
   cat("Warning: step gives NAs")
   iter=M+1
  }


 #UPDATES
 iter=iter+1
 
 funcEval=c(funcEval,f(lambda_Actual,alpha,mu))
 grad=gradient(lambda_Actual,alpha,mu,correction=correction)
 grado=c(grado, norma(grad))
 #"Distance between solutions 
 dif=norma(funcEval[length(funcEval)]-funcEval[length(funcEval)-1])
 dif_grad=norma(grado[length(grado)]-grado[length(grado)-1])
 }

   if(iter>=M)
  {
   print("Error: Number of iterations reach")
   return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-666))
  }

   #end time
   end.time <- Sys.time()
   time.taken <- abs(start.time - end.time)

   cat(rep("\n", 10)) #clean screen
   cat("\t",c("Time taken: ",time.taken," seconds"),"\n")
   #print(time.taken)
   cat("\t",c("Iterations:",iter),"\n")
   cat("\t",c("lambdas:",lambda_Actual),"\n")
   cat("\t",c("min F:",funcEval[length(funcEval)]),"\n")
   cat("\t",c("Norm(Gradient):",norma(grad)),"\n")
   cat("\t","Distance between solutions (F):",dif)
   cat("\n")
   cat("\t","Distance between gradients:",norma(grado[length(grado)]-grado[length(grado)-1]),"\n")
   if(iter<M) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=grad,error=1,dif_F=dif))#great
  
}
###########################################################################
###########################################################################

#BB_modified: BB optimization method modified

#INPUT:
#lambda: initial values of lambda
#alpha :
# mu: moments
#namefun: objetive function file name
#nameGrad: gradient function file name  
#M: number of iterations
#tolerance: tolerance of the stop criteria
#step: 4 steps to choose, number (2) is the better.
#correction: for flatter functions, this makes more convex the objetive function, 0 no correction
#1 correction
#GLL: nonmonotone line search technique (Grippo, Lampariello, and Lucidi 1986), GLL=1, or GLL=0
#for Cruz et. al. (2006) (derivative free technique)
#OUTPUT:                                                        #
#lambda_Actual=optimum lambda                                   #
# or errors: -999 (NaN values)                                  #
#            -666 (number of iterations reach)                  #
#             000 (infinite loop)                               #
#            -998 (lambdas out of range)                        #
#################################################################


BB_modified<-function(lambda, alpha, mu, nameFun = "fun.R", nameGrad = "gradient.R", M = 1000, tolerance = 1e-5,
                      step = 1, correction = 0, GLL = 1, sigma1 = 0.1, sigma2 = 0.5, phi_max = 1e20,
                      phi_min = 1e-20, N = 10)
{
 start.time <- Sys.time() #starting to measuring time
 options(digits=8)
 
 #number of points in the integration (global variable), shorting this value the convergence can be more fas
 #----------------------------------------------------------------------------------------------
 # original
 global<<- 0.0001 #0.0055 #0.001
# global<<-0.0001 #0.0055 #0.0001 #number of points in the integration (global variable)
               #shorting this value the convergence can be more fast.
 #----------------------------------------------------------------------------------------------
 #functions to use:
 #call the functions
 gradient <-source(nameGrad)$value
 f        <-source(nameFun)$value

 #Initiate variables:
 lambda_Actual=lambda
 grad=gradient(lambda,alpha,mu,correction)
 if (length(grad)==1) return(list(lambda_Opt=lambda_Actual,iter=0,grad=norma(grad),error=-998))
 grado=norma(grad)
 iter=1
 funcEval=f(lambda,alpha,mu)
 gamma=1e-4
 epsilon=1e-10
 beta=min(1,1/norma(grad))
 betas=NULL 
 m=19 # for the 4 step
 cen_max=100
 set.seed(100)
 random14=runif(cen_max)
 
 for (i in 1:length(grad))  
   if (is.na(grad[i])==T) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-999))

 dif=1
 dif_grad=1
 #STEP 1:
 while(norma(grad)>tolerance & iter<M)
 #while(dif_grad>tolerance & iter<M) #stop criteria
 {
  #STEP 2:
  if (beta<=epsilon | beta>=1/epsilon) beta=delta(grad)

  #STEP 3:
  phi=1/beta
  if (beta<=0) phi=phi_max
  if (beta>0)  phi=min(phi_max,max(phi_min,phi))
 
  #STEP 4: Nonmonotone line search to find the right lambda
  #sk=-grad search direction
  #phi= stepsize, depends of the method used 
  
  
  lambda_aux=lambda_Actual-phi*grad #iteration formula
  
  fE=f(lambda_aux,alpha,mu)

  cen=0
  while(is.na(fE)==T)
  {
   if (cen==cen_max) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-999))
   cen=cen+1
   #"Error: lambdas out of range"
   sigma=random14[cen]*(sigma2-sigma1)+sigma1
   phi=sigma*phi
   lambda_aux=lambda_Actual-phi*grad
   fE=f(lambda_aux,alpha,mu)
  } #the idea with this loop is control lambda

  if (fE==-999) 
  {
    #print("Error: lambda out of the range")
    return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=000))
    #negative values (funMEM), infinite loop
  } 

  #descending method
  index=iter-(0:min(iter-1,N))
  
  #nonmonotone line search technique (Grippo, Lampariello, and Lucidi 1986), #GLOBAL LINE SEARCH (GLL)
  if (GLL==1)
    maximo=max(funcEval[index])-gamma*(grad%*%grad)*(phi) #*(-phi) 
  #Cruz et. al. (2006) #DF-SANE (derivative free SANE)
  if (GLL==0)
  {
   nnk=grado[1]/(iter+1)^2
   maximo=max(funcEval[index])+nnk-gamma*funcEval[iter]*(phi)^2
  }
  cen=0
  while(fE>maximo) #descend condition
  {
   #STEP 5:
   cen=cen+1
   #set.seed(100)
   sigma=random14[cen]*(sigma2-sigma1)+sigma1
   phi=sigma*phi
   lambda_aux=lambda_Actual-phi*grad
   fE=f(lambda_aux,alpha,mu)
   if (cen==cen_max) 
   {
    fE=maximo #return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=000))
    cat("No disminuye", "\n")
   }
   if (is.na(fE)==T) return(list(lambda_Opt=0,iter=iter,grad=norma(grad),error=-999))
  }

  #STEP 6: update beta
  lambda_old=lambda_Actual
  lambda_Actual=lambda_aux
  gradk=gradient(lambda_Actual,alpha,mu,correction)
  if (length(gradk)==1) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-998))

  yk=gradk-grad

  #the step is controlled by this...
  if (step>4 | step<1)
  {
    cat("Error: Not valid step","\n")
    return(list(lambda_Opt=0,iter=iter,grad=norma(grad),error=999.9))
  }
  if (step==1)
  beta=-grad%*%yk/(phi*grad%*%grad)
  if (step==2)
  {
   sk=-grad
   beta=yk%*%yk/(phi*sk%*%yk)
  } 
  if (step==3)
  {
   sk=-grad
   beta2=yk%*%yk/(phi*sk%*%yk)
   beta1=-grad%*%yk/(phi*grad%*%grad)
   beta=beta1
   if (is.na(beta2)==TRUE)
   {
    cat("Warning: step gives NAs","\n")
    iter=M+1
    return(list(lambda_Opt=0,iter=iter,grad=norma(grad),error=-666))
   }
   if(beta1/beta2<=0.25) beta=beta2
  } 
  if (step==4)
  {
   sk=-grad
   beta2=yk%*%yk/(phi*sk%*%yk)
   betas[iter]=1/beta2
   sl=seq(max(1,iter-m),iter)
   beta1=min(betas[sl])
   beta=1/beta1
  } 
  if (is.na(beta)==TRUE)
  {
   cat("Warning: step gives NAs")
   iter=M+1
  }


 #UPDATES
 iter=iter+1
        if((iter %% 150)==0)
        {
             cat(c("Iteration: ", iter),"\n")
             cat(c("Min gradient so far: ", min(c(grado, norma(grad)))),"\n")
        #    cat(c("Size step: ",c,b),"\n")
        }
 
 funcEval=c(funcEval,f(lambda_Actual,alpha,mu))
 grad=gradient(lambda_Actual,alpha,mu,correction)
 if (length(grad)==1) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-998))
 grado=c(grado, norma(grad))
 #"Distance between solutions 
 dif=norma(funcEval[length(funcEval)]-funcEval[length(funcEval)-1])
 dif_grad=norma(grado[length(grado)]-grado[length(grado)-1])
 }

   if(iter>=M)
  {
   print("Error: Number of iterations reach")
   cat("Minimo gradiente de todas las iteraciones", min(grado),"\n")
   return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-666))
  }

   #end time
   end.time <- Sys.time()
   time.taken <- abs(start.time - end.time)

   cat(rep("\n", 2)) 
   cat("\t") #,c("Time taken: ",time.taken," seconds"),"\n")
   print(end.time-start.time)
   #cat(rep("\n", 1)) 
   #print(time.taken)
   cat("\t",c("Iterations:",iter),"\n")
   cat("\t",c("lambdas:",lambda_Actual),"\n")
   cat("\t",c("min F:",funcEval[length(funcEval)]),"\n")
   cat("\t",c("Norm(Gradient):",norma(grad)),"\n")
   cat("\t","Distance between solutions (F):",dif)
   cat("\n")
   cat("\t","Distance between gradients:",norma(grado[length(grado)]-grado[length(grado)-1]),"\n")
   if(iter<M) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=1,dif_F=dif,time=time.taken))#great
  
}


############################################################################

####################
#norm of a vector  #
####################

norma<-function(vector)
{
 return(sqrt(sum(vector^2)))
}


###########################################################################

#################################################################
#DESCRIPTION:                                                   #
#Value of the delta variable                                    #
#################################################################

delta<-function(grad)
{
  #norma<-source("norma.R")$value
  if (norma(grad)>1) delta=1
  if (norma(grad)>=1e-5 & norma(grad)<=1) delta=1/norma(grad)
  if (norma(grad)<1e-5) delta=1e5
 
  return(delta)
}

##############################################################################

#################################################################
#DESCRIPTION:                                                   #
#message error                                                  #
#Parameters:                                                    #
#INPUT                                                          #
#h: scalar or vector of value                                   #
#OUTPUT:                                                        #
#message                                                        #
#################################################################

error<-function(h)
{
   if((length(h)==1))
  {
   if (h==-999)
   {
    cat(rep("\n", 2)) #clean screen
    cat("+++++++++++++++++++++++++++++++++++++","\n")
    cat("+++++++++++++++++++++++++++++++++++++","\n")
    cat("Error: NaN values","\n")
   }
   if (h==-998)
   {
    cat(rep("\n", 2)) #clean screen
    cat("+++++++++++++++++++++++++++++++++++++","\n")
    cat("+++++++++++++++++++++++++++++++++++++","\n")
    cat("Error: out of range lambda values","\n")
   }
    if (h==-666)
   {
    cat(rep("\n", 2)) #borrar pantalla
    cat("+++++++++++++++++++++++++++++++++++++","\n")
    cat("+++++++++++++++++++++++++++++++++++++","\n")
    cat("Error: Number of iterations reach","\n")
   }
    if (h==0)
   {
    cat(rep("\n", 2)) #borrar pantalla
    cat("+++++++++++++++++++++++++++++++++++++","\n")
    cat("+++++++++++++++++++++++++++++++++++++","\n")
    cat("Error: Infinite loop","\n")
   }
  }
  else
  {
   #cat(rep("\n", 2)) #borrar pantalla
   cat("+++++++++++++++++++++++++++++++++++++","\n")
   cat("+++++++++++++++++++++++++++++++++++++","\n")
   cat("CONVERGENCE","\n")
  }

  return(0)
}
##############################################################################

####################################################################
#Rosembrock                                                        #
# parinit: initial guess value for the parameter to be optimized   #
# fun: Minimization function                                       #
# k: Dimensional vector size                                       #
# MaxIter: Maximum number of iterations                            #
####################################################################

Rosenbrock<-function(parinit,fun,k,MaxIter,Error,alpha,mu)
{

    # Initial point 
    lambda <- parinit 
    
    # Initial direction search
    zeta <- diag(rep(1,k)); zetaT=zeta;
    
    # Step increase in case to get a successful step(a) or one unsuccessful(b)
    a <- 3; b <- 0.5; c <- 1;
    
    # Iterativity search process until ||alfa'*zeta||<error
    error <- Error
    
    # Bucle definitions
    iter<-1; u <- matrix(1,k,k)
    points <- matrix(0,k,MaxIter)
    out <- rep(0,MaxIter)
    points[,iter] <- lambda
    out[iter] <- do.call(fun,list(lambda,alpha,mu))
    
     
    dif=1   
    while((iter<MaxIter) && (norm(u[1,],"2") > error))
    #while((iter<MaxIter) && dif > error)
    {
        vc <- 2*rep(1,k); vdir=c*rep(1,k);
        while(sum(vc)>0){
            for(i in 1:k){
                if (vc[i]>0){
                    if (do.call(fun,list((lambda+vdir[i]*zeta[,i]),alpha,mu)) <= do.call(fun,list((lambda),alpha,mu))){
                        lambda  <- lambda+vdir[i]*zeta[,i]
                        vdir[i] <- a*vdir[i]
                    }else{
                        vdir[i] <- -b*vdir[i]
                        vc[i]   <- vc[i]-1
                    }
                }
            }
        }
        iter <- iter+1
        
        points[,iter] <- lambda
        out[iter]     <- do.call(fun,list((lambda),alpha,mu)) 
        dif=c(dif,abs(out[iter]-out[iter-1]))  
        # Computation of new search directions
        vdif <- points[,iter]-points[,(iter-1)]
        if(sum(vdif^2)>0){
            u <- diag(rep(0,k))
            for(i in 1:k){
                for(j in i:k){
                    u[i,] <- u[i,]+vdif[j]*zetaT[,j]
                }
            }
            for(i in 1:k){
                w <- u[i,]
                if(i>1){
                	for(j in 1:(i-1)){
                    	w <- w-(u[i,]*zeta[,j])*zeta[,j] 
                	}
                }
                if(norm(w,"2")==0){
                    zeta <- zetaT
                    break
                }else{
                    zeta[,i] <- w/norm(w,"2")
                }
            }
            zetaT <- zeta
        }else{
            a <- a/2; b <- b/2
            if(a<1){
               a <- 1 
            }
            c <- c/2
        }
        if((iter %% 50)==0){
             cat(c("Iteration: ", iter),"\n")
        #    cat(c("Min value so far: ", out[iter]),"\n")
        #    cat(c("Size step: ",c,b),"\n")
        }
    }
    VecMin <- points[,iter]
    OutMin <- do.call(fun,list(VecMin,alpha,mu))
    ArrVec <- points[,1:iter]
    ArrOut <- out[1:iter]
    
    cat(rep("\n", 10)) #clean screen
    cat(c("Process finished at iteration: ", iter),"\n") 
    cat(c("Min value found: ", out[iter]),"\n")
    cat(c("Differences between solutions: ", dif[iter]),"\n")
    cat(c("Minimum lambdas: ", VecMin),"\n")
 

 return(list(VecMin=VecMin,OutMin=OutMin,ArrVec=ArrVec,ArrOut=ArrOut))
}

#############################################################################

####################################################################
#maxent_parameters:                                                #
#Calculates the input paramentes of SME using the simulated data   #
#                                                                  #
#INPUT                                                             #
# file: the name of the file with the simulated data               #
# k: Dimensional vector size                                       #
#OUTPUT                                                            #
# alpha                                                            #
# moments                                                          #
####################################################################

alpha_function <- function(k){
  dec   <- 1.5
  dec1  <- rep(dec,k)
  kk1   <- 1:k
  alpha <- dec1/kk1
  return(alpha)
}

maxent_parameters <- function(S, k)
{

 M=length(S)
 #calculation of alpha
 alpha <- alpha_function(k)

 #calculation of momnets mu
 
 #Prob of N=0
 #po=length(n[n==0])/M
 po=length(S[S==0])/M

 laplace=rep(0,k)
 mu=rep(0,k)

 #laplace transform
 for (i in 1:k)
  laplace[i]=(1/M)*sum(exp(-alpha[i]*S)) 

 #moments
 mu=(laplace-po)/(1-po)

 return(list(mu=mu))

}
##########################################################################

MEM_parameters<-function(S,k,N,lambda_initial)
{
 M=length(S)
 if (length(N)==length(lambda_initial)) N=N+1 #se pudiera quitar
 #calculation of alpha
 dec=1.5
 dec1=rep(dec,k)
 kk1=1:k
 alpha_no0=dec1/kk1
 alpha=c(0,alpha_no0)

 mA=source("matrixA.R")$value
 A1=mA(alpha_no0,N)
 A=rbind(rep(1,N),A1)
 #calculation of moments mu
 
 #Prob of N=0
 #po=length(n[n==0])/M
 po=length(S[S==0])/M

 laplace=rep(0,k)
 mu=rep(0,k)

 #laplace transform
 for (i in 1:k)
 laplace[i]=(1/M)*sum(exp(-alpha_no0[i]*S)) 

 #moments
 mu=(laplace-po)/(1-po)
 mu=c(1,mu)


 #lambda
 lambda0=log(sum(exp(-t(A1)%*%lambda_initial))) #*OJO: caso Poisson, necesita un (-)
 lambda=c(lambda0,lambda_initial)

 
 return(list(mu=mu,alpha=alpha,lambda=lambda))

}


##########################################################################
##########################################################################

bootstrap_mu<-function(k,S,rep=1000,seed=10,p1=0.20,p2=0.25)
{
 M=length(S)
 set.seed(seed)

 #calculation of alpha
 dec=1.5
 dec1=rep(dec,k)
 kk1=1:k
 alpha=dec1/kk1
 
 po=length(S[S==0])/length(S)
 
 laplace=NULL
 mu_f=NULL
 #laplace transform
 for (i in 1:k)
  laplace[i]=(1/M)*sum(exp(-alpha[i]*S)) 
 #moments
 mu_f=(laplace-po)/(1-po) 

 #####################################################
  mu_matrix=matrix(0,ncol=k,nrow=rep)
  for (j in 1:rep)
  {
    size=trunc(runif(1,p1*M,p2*M)) #muestra entre 10% y 15%
    S_sample=NULL #INICIALIZAR
    S_sample<-sample(S, size, replace = FALSE, prob = NULL)

    #calculation of moments of mu
 
    #Prob of N=0
    po=length(S_sample[S_sample==0])/length(S_sample)

    laplace=rep(0,k)
    mu=rep(0,k)

    #laplace transform
    for (i in 1:k)
      laplace[i]=(1/length(S_sample))*sum(exp(-alpha[i]*S_sample)) 

    #moments
    mu=(laplace-po)/(1-po)
    mu_matrix[j,1:k]=mu[1:k]
   }

   #interval at 95% of confidence
   mu_sup=NULL
   mu_inf=NULL
   vect=NULL
   for (i in 1:k)
   {
    vect=sort(mu_matrix[,i])
    #percentiles del 97.5%
    pos_97<-round(rep*(0.975))
    #percentiles del 2.5%
    pos_2<-round(rep*(0.025))
   
    mu_sup[i]=vect[pos_97]
    mu_inf[i]=vect[pos_2]
   }

  int_mu=cbind(mu_inf,mu_sup)

  return(list(alpha=alpha,int_mu=int_mu,mu=mu_f))
}
###########################################################################
###########################################################################
#rep: number of repetitions
#p1 & p2: value for sellection of the sample
bootstrap_mu_gamma<-function(k,S,rep=1000,seed=10,p1=1,p2=1)
{
 M=length(S)
 set.seed(seed)

 #calculation of alpha
 dec=1.5
 dec1=rep(dec,k)
 kk1=1:k
 alpha=dec1/kk1
 

 #MOMENTS OF THE DATA
 po=length(S[S==0])/length(S) 
 laplace=NULL
 mu_f=NULL
 #laplace transform
 for (i in 1:k)
  laplace[i]=(1/M)*sum(exp(-alpha[i]*S)) 
 #moments
 mu_empirical=(laplace-po)/(1-po) 

 #THEORETICAL MOMENTS
 mu1<-source("moments_gamma.R")$value
 mu_theorical=mu1(alpha)


 #####################################################
  mu_matrix=matrix(0,ncol=k,nrow=rep)
  for (j in 1:rep)
  {
    size=trunc(runif(1,p1*M,p2*M)) #p1=1,   p2=1 takes all the sample
                                   #p1=0.5, p2=1 takes between 50% and 100%
    S_sample=NULL #INICIALIZAR
    S_sample<-sample(S, size, replace = TRUE, prob = NULL)

    #calculation of moments of mu
 
    #Prob of N=0
    po=length(S_sample[S_sample==0])/length(S_sample)

    laplace=rep(0,k)
    mu=rep(0,k)

    #laplace transform
    for (i in 1:k)
      laplace[i]=(1/length(S_sample))*sum(exp(-alpha[i]*S_sample)) 

    #moments
    mu=(laplace-po)/(1-po)
    mu_matrix[j,1:k]=mu[1:k]
   }

   #interval at 95% of confidence
   mu_sup=NULL
   mu_inf=NULL
   vect=NULL
   for (i in 1:k)
   {
    vect=sort(mu_matrix[,i])
    #percentiles del 97.5%
    pos_97<-round(rep*(0.975))
    #percentiles del 2.5%
    pos_2<-round(rep*(0.025))
   
    mu_sup[i]=vect[pos_97]
    mu_inf[i]=vect[pos_2]
   }

  int_mu=cbind(mu_inf,mu_sup)

  return(list(alpha=alpha,int_mu=int_mu,mu_empirical=mu_empirical,
  mu_theorical=mu_theorical))
}
############################################################################

theorical_interval_mu_gamma<-function(ele, a, b, M, k,rep=1000,seed=10)
{
 set.seed(seed)
 #calculation of alpha
 dec=1.5
 dec1=rep(dec,k)
 kk1=1:k
 alpha=dec1/kk1


 mu_matrix=matrix(0,ncol=k,nrow=rep)
 for (j in 1:rep)
 {
  #dist Poisson
  n=rpois(M,ele) #N

  S=vector(mode = "numeric", length = M)
  # for each k 
  for (kk in 1:M)
  {
   xc=(rgamma(n[kk], shape = a, scale = b)) 
   S[kk]=sum(xc)
  }

    hist(S)
    #Prob of N=0
    po=length(S[S==0])/length(S)
    laplace=rep(0,k)
    mu=rep(0,k)
    #laplace transform
    for (i in 1:k)
      laplace[i]=(1/length(S))*sum(exp(-alpha[i]*S)) 
    #moments
    mu=(laplace-po)/(1-po)

   mu_matrix[j,]=mu
  }

   means=NULL
   sds=NULL
    for (i in 1:k)
   {
    means[i]=mean(mu_matrix[,i])
    sds[i]=sd(mu_matrix[,i])
   }

  mu_inf=means-qnorm(1-(0.05)/2)*(sds/sqrt(rep))
  mu_sup=means+qnorm(1-(0.05)/2)*(sds/sqrt(rep))
  int_mu=cbind(mu_inf,mu_sup)

  return(list(alpha=alpha,int_mu=int_mu))
}
 
############################################################################

theorical_interval_mu_simple_log<-function(alpha, me, sd, M, k, gamma=.95,
                                           seed=1)
{
    set.seed(1)
    S=rlnorm(M, meanlog = me, sdlog = sd)
    hist(S)
    #Prob of N=0
    po=length(S[S==0])/length(S)
    laplace=rep(0,k)
    laplace2=rep(0,k)
    mu=rep(0,k)
    #laplace transform
    for (i in 1:k)
    {
      h=exp(-alpha[i]*S)
      laplace[i]=(1/length(S))*sum(h) 
      laplace2[i]=(1/length(S))*sum(h^2) 
    }
    #moments
    mu=(laplace-po)/(1-po)
    sds=sqrt(laplace2-(laplace)^2)
    print(sds)

  mu_inf=mu-qnorm(1-(1-gamma)/2)*(sds/sqrt(M))
  mu_sup=mu+qnorm(1-(1-gamma)/2)*(sds/sqrt(M))
  #mu_inf=mu-0.1*(sds/sqrt(M))
  #mu_sup=mu+0.1*(sds/sqrt(M))
  int_mu=cbind(mu_inf,mu_sup)

  return(list(int_mu=int_mu))
}

###########################################################################


#CONJUGATE GRADIENT METHOD

BB2<-function(lambda,alpha,mu,nameFun="fun.R",nameGrad="gradient.R",M=1000,tolerance=1e-6,step=1,
              correction=0,phi_max=1e20,phi_min=1e-20)
{
 options(digits=8)
 #functions to use:
 #call the functions
 gradient<-source(nameGrad)$value
 f<-source(nameFun)$value
 #norma<-source("norma.r")$value
 #delta<-source("delta.r")$value

 #Initiate variables:
 lambda_Actual=lambda
 grad=gradient(lambda,alpha,mu,correction)
 grado=norma(grad)
 iter=1
 funcEval=f(lambda,alpha,mu)
 gamma=1e-4
 epsilon=1e-10
 sigma1=0.1
 sigma2=0.5
 #sigma1=1e-4
 #sigma2=0.9
 beta=min(1,1/norma(grad))
 N=10
 m=19 # for the 4 step
 set.seed(1)
 random14=runif(100)
 
 for (i in 1:length(grad))  
 if (is.na(grad[i])==T) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-999))

 dk_aux=NULL
 dk=NULL
 dif=1
 dif_grad=1
 #STEP 1:
 while(norma(grad)>tolerance & iter<M)
 {
  #STEP 2:
  if (is.na(beta)==TRUE)  
  {
   cat("Beta INF","\n")
   beta=delta(grad)
   #return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-999))
  }

  if (beta<=epsilon | beta>=1/epsilon) beta=delta(grad)

  #STEP 3:
  phi=1/beta
  if (beta<=0) phi=phi_max
  if (beta>0)  phi=min(phi_max,max(phi_min,phi))


  ###
  if (iter==1)
     dk=-grad
  if (iter>1)
  {
     dk_aux=dk
     Bk=(norma(grad)^2)/(abs(grad%*%dk_aux)+(grado[iter-1])^2)
     #Bk=max(c(0,min(c((grad%*%yk)/(grado[iter-1])^2),(norma(grad)^2)/(grado[iter-1])^2)))    
     dk=-grad+Bk*dk_aux    
  }
  
  #STEP 4: Nonmonotone line search to find the right lambda
  #dk=-grad search direction
  #phi= stepsize, depends of the method used 
  
  lambda_aux=lambda_Actual+phi*dk #iteration formula
  
  fE=f(lambda_aux,alpha,mu)

  cen=0
  while(is.na(fE)==T)
  {
   if (cen==100) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-999))
   cen=cen+1
   #"Error: lambdas out of range"
   sigma=random14[cen]*(sigma2-sigma1)+sigma1
   phi=sigma*phi

   lambda_aux=lambda_Actual+phi*dk #**
   fE=f(lambda_aux,alpha,mu)
  } #the idea with this loop is control lambda

  if (fE==-999) 
  {
    #print("Error: lambda out of the range")
    return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=000))
    #negative values (funMEM), infinite loop
  } 


  #descending method
  index=iter-(0:min(iter-1,N))
  maximo=max(funcEval[index])-gamma*(grad%*%grad) #*phi
  cen=0
  while(fE>maximo) #descend condition
  {
   #STEP 5:
   cen=cen+1
   sigma=random14[cen]*(sigma2-sigma1)+sigma1
   phi=sigma*phi
   lambda_aux=lambda_Actual+phi*dk #**
   fE=f(lambda_aux,alpha,mu)
   if (is.na(fE)==T) return(list(lambda_Opt=0,iter=iter,grad=norma(grad),error=-999))
   if (cen==100)fE=maximo # return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=000))
  }

  #STEP 6: update beta
  lambda_old=lambda_Actual
  lambda_Actual=lambda_aux
  yk=gradient(lambda_Actual,alpha,mu,correction)-grad
  #print(yk)

  #the step is controlled by this...
  if (step>2 | step<1)
  {
    cat("Error: Not valid step","\n")
    return(list(lambda_Opt=0,iter=iter,grad=norma(grad),error=999.9))
  }
  if (step==1)
  beta=-grad%*%yk/(phi*grad%*%grad)
  if (step==2)
  {
   sk=-grad
   beta=yk%*%yk/(phi*sk%*%yk)
  } 

 #UPDATES
 iter=iter+1
 if((iter %% 150)==0)
 {
  cat(c("Iteration: ", iter),"\n")
  cat(c("Min gradient so far: ", norma(grad)),"\n")
  #cat(c("Size step: ",c,b),"\n")
 }
 
 funcEval=c(funcEval,f(lambda_Actual,alpha,mu))
 grad=gradient(lambda_Actual,alpha,mu,correction)
 grado=c(grado, norma(grad))
 #"Distance between solutions 
 dif=norma(funcEval[length(funcEval)]-funcEval[length(funcEval)-1])
 dif_grad=norma(grado[length(grado)]-grado[length(grado)-1])
 }

   if(iter==M)
  {
   #print("Error: Number of iterations reach")
   return(list(lambda_Opt=lambda_Actual,iter=iter,grad=norma(grad),error=-666))
  }

  
   cat(rep("\n", 10)) #clean screen
   cat("\t",c("Iterations:",iter),"\n")
   cat("\t",c("lambdas:",lambda_Actual),"\n")
   cat("\t",c("min F:",funcEval[length(funcEval)]),"\n")
   cat("\t",c("Norm(Gradient):",norma(grad)),"\n")
   cat("\t","Distance between solutions (F):",dif)
   cat("\n")
   cat("\t","Distance between gradients:",norma(grado[length(grado)]-grado[length(grado)-1]),"\n")
   if(iter<M) return(list(lambda_Opt=lambda_Actual,iter=iter,grad=grad,error=1,dif_F=dif))#great
  
}




