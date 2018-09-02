#sample_size_variability_functions: the functions used to generate many times the denisties and calculate sampling variability!!!!!! 

######################################################################################################
simulR<-function(S,tam=50,Mnumber=1000,kMoments=8,Tol=1e-4,lambda0=rep(1,kMoments))
{
  if (length(lambda0)!=kMoments)
  {
    print("Error")
    return(0)
  }
  resultados=NULL
  for (i in 1:Mnumber)
  {
    graphics.off()
    #vector i
    set.seed(i)
    indS=NULL
    indS=sample(1:length(S),size=tam)
    S9=NULL
    S9=S[indS]
    
    cat(c("**Iteration: ", i),"\n") #indicate where I am.
    
    
    k=kMoments
    mp=maxent_parameters(S9,k) #aqui se deben incluir los ceros (OJO)
    mu=mp$mu
    alpha=mp$alpha
    
    global<<-0.0001 #0.0055 #0.0001 #number of points in the integration (global variable)
    #shorting this value the convergence can be more fast.
    
    #fix always a tolerance which gives good results in general, in some cases
    #you can improved, but not always. Depends of the case.
    tol=Tol
    lambda=lambda0
    source("opt_methods.r")
    #Initial values
    h=BB_modified(lambda,alpha,mu,M=10000,tolerance=tol,step=2,correction=0,GLL=1,
                  sigma1=0.1,sigma2=0.5,phi_max=1e30,phi_min=1e-30,N=100)
    if (h$error==000 | h$error==-999 | h$error==-666) error(h$error)
    lambda=h$lambda_Opt
    
    resultados=rbind(resultados,lambda)
  }
  
  return(resultados)
}