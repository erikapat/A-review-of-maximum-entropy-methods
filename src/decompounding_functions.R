
#-----------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

###########################################################################

#GENERAL: FOR POISSON, Geometric, Negative Binomial
#gives the laplace transformation of the individual losses
laplace_individuals <- function(lambda=NULL,alpha=NULL,distribution="Poisson",...,globalZ=1e-3,mu=NULL)
{
  k=length(lambda)
  par<-list(...)
  if (distribution!="geometric" & distribution!="Poisson" & distribution!="negative binomial" &  distribution!="Binomial") {
    #if (length(distribution)==0)
    cat("ERROR: frequency distribution not in the (a,b,0) family",'\n')
    return(0)
  } 
  
  vals=c(length(lambda),length(alpha),length(mu))
  if (all(vals==0)){
    cat("ERROR: NULL parameters, please check",'\n')
    return(0)
  }
  
  #laplace function of S
  laplace_S=NULL
  if (length(mu)==0){
    int=NULL
    y<-seq(from=0,to=1,by=globalZ)
    f_y=sum_exp(lambda,alpha,y)
    Z<-sintegral(y,f_y)$value #Z normalization parameter
    #Z<-trapz(y,f_y)
    f_asterico=f_y/Z
    for (i in 1:k)
      int[i]<-sintegral(y,(y^alpha[i])*f_asterico)$value #Z normalization parameter
  }
  
  #laplace of the aggregated
  if (distribution=="geometric")
  {
    distribution="negative binomial"
    par$size=1
  }
  
  if (distribution=="Poisson")
  {
    if (!"Lambda" %in% names(par))
      stop("value of 'Lambda' missing")
    
    l=par$Lambda
    
    po=exp(-l)
    if (length(lambda)==0 & length(alpha)==0)
      laplace_S=mu
    if (length(mu)==0)
      laplace_S=po+(1-po)*int
    
    
    #laplace individual losses
    laplace_l=(1/(l))*log(laplace_S)+1
    
    
  }
  if (distribution=="negative binomial")
  {
    if (!"prob" %in% names(par))
    {
      
      cat("value of 'prob' missing",'\n')
      return(0)
    }
    if (!"size" %in% names(par))
    {
      cat("value of 'size' missing",'\n')
      return(0)
    }
    beta=1/par$prob - 1
    r= par$size
    
    #laplace of the aggregated
    po=(1/(1+beta))^r
    if (length(lambda)==0 & length(alpha)==0)
      laplace_S=mu
    if (length(mu)==0)
      laplace_S=po+(1-po)*int
    
    laplace_l=(1-laplace_S^(-1/r))/beta + 1 #BIN NEGATIVE
    
  }
  if (distribution=="Binomial")
  {
    if (!"prob" %in% names(par))
    {
      
      cat("value of 'prob' missing",'\n')
      return(0)
    }
    if (!"size" %in% names(par))
    {
      cat("value of 'size' missing",'\n')
      return(0)
    }
    po=par$prob
    n= par$size
    
    #laplace of the aggregated
    
    if (length(lambda)==0 & length(alpha)==0)
      laplace_S=mu
    if (length(mu)==0)
      laplace_S=po+(1-po)*int
    
    laplace_l=(laplace_S^(1/n)-(1-po))/po #BIN NEGATIVE
    
  }
  
  hu=laplace_l
  
  return(hu)
}

