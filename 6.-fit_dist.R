

#ProcesoDeSimulacionAutomation.R
#SIMULATIONS

#rm(list=ls()) #remover cualquier variable del espacio de trabajo :)
#setwd("C:\\Users\\Erika\\Dropbox\\4ErikaIII (continuacion)\\paper (5) scaling\\code")

main <- function(TAM = 10, minVal=3.5,maxVal=7,interVals=30) {
  
  library(actuar)
  library(fitdistrplus)
  #sizes
  M1 <- TAM
  
  tail.densi <- NULL
  cov.densi <- NULL
  cov.cdf  <- NULL
  tail.cdf <- NULL
  for (i in 1:200){
    set.seed(i) 
    cat('[INFO]:', ' Iteration', i, '\n')
    # without correlations
    n1=rpois(M1, 80)
    n2=rpois(M1, 60)
    n3=rbinom(M1, 70, 0.5) #mean 19.5
    n4=rbinom(M1, 62, 0.5)  #mean 27
    n5=rbinom(M1, 50, 0.5)   #mean 7.5
    n6=rbinom(M1, 76, 0.5)  #mean 32
    n7=rnbinom(M1,80, 0.3)
    if (M1 < 500)
      n8=rnbinom(1000, 90, 0.8) 
    else
      n8=rnbinom(M1, 90, 0.8) 
    
    alpha=20
    M=85
    c=15
    source("champernowne.R")
    S1=vector(mode = "numeric", length = M1)
    indi=NULL
    # for each k 
    jk1=NULL
    for (k in 1:M1)
    {
      xc=NULL
      xc=champer_random(n1[k],alpha,c,M) 
      S1[k]=sum(xc)
      jk1=c(jk1,xc)
    }
    hist(S1,breaks=25,freq=FALSE,include.lowest = TRUE)
    
    
    
    #-------------------------------------------------------------------------------------
    x <- seq(min(S1), max(S1), length= 100)
    f1 <- fitdist(S1[S1>0],"norm")
    densi.f1 <- dnorm(x, f1$estimate[1], f1$estimate[2])
    lines(x,densi.f1, col ='green')
    
    cdfcomp(list(f1),
            xlogscale = TRUE, ylogscale = TRUE,
            legendtext = c("norm"))
    
    #---------------------------------------------------------------------------------------
    
    library(evd)
    mn=-0.01
    sdl=2
    S2=vector(mode = "numeric", length = M1)
    indi=NULL
    # for each k 
    jk1=NULL
    for (k in 1:M1)
    {
      xc=NULL
      xc=(rlnorm(n2[k], meanlog = mn, sdlog = sdl)) 
      S2[k]=sum(xc)
    }
    hist(S2,breaks=25,freq=FALSE,include.lowest = TRUE)
    min(S2)
    max(S2)
    
    #--------------------------------------------------------------------------------------------
    x <- seq(min(S2), max(S2), length= 100)
    f2 <- fitdist(S2[S2>0],"lnorm")
    densi.f2 <- dlnorm(x, f2$estimate[1], f2$estimate[2])
    lines(x, densi.f2, col ='green')
    
    cdfcomp(list(f2),
            xlogscale = TRUE, ylogscale = TRUE,
            legendtext = c("lnorm"))
    
    #-------------------------------------------------------------------------------------------
    
    
    mn3=10
    sdl3=85
    
    S3=vector(mode = "numeric", length = M1)
    indi=NULL
    # for each k 
    for (k in 1:M1)
    {
      xc=NULL
      xc=rpareto(n3[k], mn3, scale=sdl3)
      S3[k]=sum(xc)
    }
    #GRAPH of the simulation
    hist(S3[S3>0],breaks=15,freq=FALSE,include.lowest = TRUE)
    
    #--------------------------------------------------------------------------------------------
    x <- seq(min(S3), max(S3), length= 100)
    f3 <- fitdist(S3[S3>0], "gamma")
    densi.f3 <- dgamma(x, f3$estimate[1], f3$estimate[2])
    lines(x, densi.f3, col ='green')
    
    cdfcomp(list(f3),
            xlogscale = TRUE, ylogscale = TRUE,
            legendtext = c("gamma"))
    
    #-------------------------------------------------------------------------------------------
    
    alpha=10
    M=125
    c=45
    source("champernowne.R")
    S4=vector(mode = "numeric", length = M1)
    indi=NULL
    # for each k 
    for (k in 1:M1)
    {
      xc=NULL
      xc=champer_random(n4[k],alpha,c,M) 
      S4[k]=sum(xc)
    }
    hist(S4,breaks=25,freq=FALSE,include.lowest = TRUE)
    
    #--------------------------------------------------------------------------------------------
    x <- seq(min(S4), max(S4), length= 100)
    f4 <- fitdist(S4[S4>0], "norm")
    densi.f4 <- dnorm(x, f1$estimate[1], f1$estimate[2])
    lines(x, densi.f4, col ='green')
    
    cdfcomp(list(f4),
            xlogscale = TRUE, ylogscale = TRUE,
            legendtext = c("norm"))
    
    #-------------------------------------------------------------------------------------------
    
    
    S5=vector(mode = "numeric", length = M1)
    indi=NULL
    # for each k 
    for (k in 1:M1)
    {
      xc=NULL
      xc=rgamma(n5[k],4500, 15)
      S5[k]=sum(xc)
    }
    hist(S5,breaks=25,freq=FALSE,include.lowest = TRUE)
    
    #--------------------------------------------------------------------------------------------
    x <- seq(min(S5), max(S5), length= 100)
    f5 <- fitdist(S5[S5>0], "norm")
    densi.f5 <- dnorm(x, f5$estimate[1], f5$estimate[2])
    lines(x, densi.f5, col ='green')
    
    cdfcomp(list(f5),
            xlogscale = TRUE, ylogscale = TRUE,
            legendtext = c("norm"))
    
    #-------------------------------------------------------------------------------------------
    
    
    S6=vector(mode = "numeric", length = M1)
    indi=NULL
    # for each k 
    for (k in 1:M1)
    {
      xc=NULL
      xc=rgamma(n6[k],9000, 35)
      S6[k]=sum(xc)
    }
    hist(S6,breaks=25,freq=FALSE,include.lowest = TRUE)
    
    #--------------------------------------------------------------------------------------------
    x <- seq(min(S6), max(S6), length= 100)
    f6 <- fitdist(S6[S6>0], "norm")
    densi.f6 <- dnorm(x, f6$estimate[1], f6$estimate[2])
    lines(x, densi.f6, col ='green')
    
    cdfcomp(list(f6),
            xlogscale = TRUE, ylogscale = TRUE,
            legendtext = c("norm"))
    
    #-------------------------------------------------------------------------------------------
    
    
    S7=vector(mode = "numeric", length = M1)
    indi=NULL
    # for each k 
    for (k in 1:M1)
    {
      xc=NULL
      xc=rweibull(n7[k],200, 50)
      S7[k]=sum(xc)
    }
    hist(S7,breaks=25,freq=FALSE,include.lowest = TRUE)
    
    #--------------------------------------------------------------------------------------------
    x <- seq(min(S7), max(S7), length= 100)
    f7 <- fitdist(S7[S7>0], "lnorm")
    densi.f7 <- dlnorm(x, f7$estimate[1], f7$estimate[2])
    lines(x, densi.f7, col ='green')
    
    cdfcomp(list(f7),
            xlogscale = TRUE, ylogscale = TRUE,
            legendtext = c("lnorm"))
    
    #-------------------------------------------------------------------------------------------
    
    
    mn8=5.5
    sdl8=5550
    if (M1 < 500)
      M8=1000 
    else
      M8=M1 
    S8=vector(mode = "numeric", length = M8)
    indi=NULL
    # for each k 
    for (k in 1:M8)
    {
      xc=NULL
      xc=rpareto(n8[k], mn8, scale=sdl8)
      S8[k]=sum(xc)
    }
    hist(S8,breaks=25,freq=FALSE,include.lowest = TRUE)
    um <- 43000
    hist(S8[S8>um],breaks=15,freq=FALSE,include.lowest = TRUE)
    
    
    #--------------------------------------------------------------------------------------------
    x <- seq(min(S8[S8>um]), max(S8[S8>um]), length= 100)
    f8 <- fitdist(S8[S8>um], "lnorm")
    densi.f8 <- dlnorm(x, f8$estimate[1], f8$estimate[2])
    lines(x, densi.f8, col ='green')
    
    cdfcomp(list(f8),
            xlogscale = TRUE, ylogscale = TRUE,
            legendtext = c("lnorm"))
    
    #-------------------------------------------------------------------------------------------
    
    #GRAPH of the simulation
    
    Sh=S1+S2+S3+S4+S5+S6+S7
    S=c(Sh,S8[S8>um])	
    S=S[S>32000]
    hist(S,breaks=35,freq=FALSE,include.lowest = TRUE)
    summary(S)
    
    S=S/1e4
    hist(S,breaks=35,freq=FALSE,include.lowest = TRUE)
    summary(S)
    
    #-------------------------------------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------------------------------------
    
    # with the convolution package
    
    library(distr)
    ff1 <- Norm(f1$estimate[1], f1$estimate[2])          # N is signature for normal dist
    ff2 <- Lnorm(f2$estimate[1], f2$estimate[2]) # same for log-normal
    ff3 <- Gammad(f3$estimate[1], f3$estimate[2]) 
    ff4 <- Norm(f4$estimate[1], f4$estimate[2])
    ff5 <- Norm(f5$estimate[1], f5$estimate[2])
    ff6 <- Norm(f6$estimate[1], f6$estimate[2])
    ff7 <- Lnorm(f7$estimate[1], f7$estimate[2])
    ff8 <- Lnorm(f8$estimate[1], f8$estimate[2])
    
    conv <- convpow(ff1+ff2 + ff3+ ff4 + ff5 + ff6 + ff7,1)             # object of class AbscontDistribution
    f.Z   <- d(conv)                    # distribution function
    Ff.Z  <- p(conv)
    
    
    z <- seq(minVal*1e4, maxVal*1e4, length = interVals) 
    z <- Vectorize(z)
    Z <- S1 +S2 +S3 +S4 +S5 +S6 + S7+ S8
    hist(Z,freq=F)
    lines(z,f.Z(z),lty=2, type = 'l', col="red")
    
    
    S <- read.table(file = "S_2ndExample.dat")
    S <- S[,1]
    S <- S[S>3.7]
    hist(S,freq=F, ylim =c(0,2))
    lines(spline(z*1e-4  + .33,f.Z(z)*1e4),lty=2, type = 'l', col="red")
    lines(spline(z[z>40700]*1e-4  + .33,f.Z(z)[z>40700]*1e4),lty=2, type = 'l', col="blue")
    
    x <- seq(55000, 80000, length= 100)
    densi.f8 <- dlnorm(x, f8$estimate[1], f8$estimate[2])
    cdf.f8 <- plnorm(x, f8$estimate[1], f8$estimate[2])
    lines(spline(x*1e-4-1,densi.f8*1e4),lty=2, type = 'l', col="blue")
    
    tail.densi <- rbind(tail.densi, densi.f8*1e4)
    tail.cdf   <- rbind(tail.densi, cdf.f8*1e4)
    cov.densi  <- rbind(cov.densi, f.Z(z)*1e4)
    cov.cdf    <- rbind(cov.cdf, Ff.Z(z))
  }
  cov.densi  <- rbind(z*1e-4  + .33, cov.densi)
  tail.densi <- rbind(x*1e-4-1,      tail.densi)
  cov.cdf  <- rbind(z*1e-4  + .33,   cov.cdf)
  tail.cdf <- rbind(x*1e-4-1,        tail.cdf)
  
  ifile <- paste0("FitDistributions/convolution_density_", TAM, '.dat')
  cat('[INFO]:', ' Saving', ifile, '\n')
  write.table(cov.densi , file = ifile, quote = F, row.names = F, col.names = F)
  
  ifile <- paste0("FitDistributions/tail_density_", TAM, '.dat')
  cat('[INFO]:', ' Saving', ifile, '\n')
  write.table(tail.densi , file = ifile, quote = F, row.names = F, col.names = F)
  
  ifile <- paste0("FitDistributions/convolution_CDF_", TAM, '.dat')
  cat('[INFO]:', ' Saving', ifile, '\n')
  write.table(cov.cdf , file = ifile, quote = F, row.names = F, col.names = F)
  
  ifile <- paste0("FitDistributions/tail_CDF_", TAM, '.dat')
  cat('[INFO]:', ' Saving', ifile, '\n')
  write.table(tail.cdf , file = ifile, quote = F, row.names = F, col.names = F)
}

#---------------------------------------------------------------------------------------------------------
options(echo = TRUE)
argv <- commandArgs(trailingOnly = T)
print(argv)

if (length(argv) < 1) {
  cat("Params: \n")
  stop()
} #else if (length(argv) > 2) {
#cat("Params:  \n")
# stop()
#}

system.time(state <- main(as.integer(argv[1]), as.numeric(argv[2]), as.numeric(argv[3]), 
                          as.integer(argv[4])))
