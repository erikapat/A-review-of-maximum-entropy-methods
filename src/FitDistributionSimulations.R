
# Burr Distribution 
simulFitBurr <-function(S,tam = 50,Mnumber = 200, minVal = 3.5,maxVal = 7,interVals = 30){

  accum.dist <- NULL
  densidad <- NULL
  
  if (tam > 200){
    l.value <- 0.251
  }else{
    l.value <- 0.258
  }
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
    

    fendo.B <- fitdist(S9[S9>0], "burr", start = list(shape1 = 0.9, shape2 = 2,rate = .1),
                       lower = c(0.3, .5,  l.value)) 
    densi   <- dburr(seq(minVal, maxVal, length = interVals), shape1 = fendo.B$estimate[1], shape2 = fendo.B$estimate[2], rate = fendo.B$estimate[3],
                log = FALSE)
    cdf.B  <- pburr(seq(minVal, maxVal, length = interVals), shape1 = fendo.B$estimate[1], shape2 = fendo.B$estimate[2], rate = fendo.B$estimate[3],
                     log = FALSE)
    
    densidad   <- cbind(densidad, densi)
    accum.dist <- cbind(accum.dist, cdf.B)
  }
  x = seq(minVal, maxVal, length = interVals)
  densi.estimated <- data.frame(x = x, densidad )
  
  return(densi.estimated)
  
}

#--------------------------------------------------------------------------------------------------
# GEV

simulFitGEV <-function(S,tam=50,Mnumber=200, minVal=3.5, maxVal=7, interVals=30){
  accum.dist <- NULL
  densidad <- NULL
  for (i in 1:Mnumber){
    graphics.off()
    #vector i
    set.seed(i)
    indS=NULL
    indS=sample(1:length(S),size=tam)
    S9=NULL
    S9=S[indS]
    
    cat(c("**Iteration: ", i),"\n") #indicate where I am.
    
    par <- gev.fit(S9[S9>0], show = FALSE)
    densi  <- dgev(seq(minVal, maxVal, length = interVals), xi = par$mle[3], mu = par$mle[1], beta = par$mle[2])
    cdf.B  <- pgev(seq(minVal, maxVal, length = interVals), xi = par$mle[3], mu = par$mle[1], beta = par$mle[2])
    
    densidad   <- cbind(densidad, densi)
    accum.dist <- cbind(accum.dist, cdf.B)
  }
  x = seq(minVal, maxVal, length = interVals)
  densi.estimated <- data.frame(x = x, densidad )
  
  return(densi.estimated)
  

}
#-----------------------------------------------------------------------------------------------

# LogNormal Distribution 
simulFitLogNormal <-function(S,tam = 50,Mnumber = 200, minVal = 3.5,maxVal = 7,interVals = 30){
  
  accum.dist <- NULL
  densidad <- NULL
  for (i in 1:Mnumber)
  {
    graphics.off()
    #vector i
    set.seed(i*5)
    indS=NULL
    indS=sample(1:length(S),size=tam)
    S9=NULL
    S9=S[indS]
    
    cat(c("**Iteration: ", i),"\n") #indicate where I am.
    
    fendo.B <- fitdist(S9[S9>0], "lnorm", method = 'mle')
    densi   <- dlnorm(seq(minVal, maxVal, length = interVals), fendo.B$estimate[1], fendo.B$estimate[2])
    cdf.B  <-  plnorm(seq(minVal, maxVal, length = interVals), fendo.B$estimate[1], fendo.B$estimate[2])
    
    densidad   <- cbind(densidad, densi)
    accum.dist <- cbind(accum.dist, cdf.B)
  }
  x = seq(minVal, maxVal, length = interVals)
  densi.estimated <- data.frame(x = x, densidad )
  
  return(densi.estimated)

}

#-------------------------------------------------------------------------------------------------

# GRAPH
simGraph<-function(densi,S,minVal=0,maxVal=2,porc=1.83)
{
  
  interVals=dim(densi)[2]
  x=seq(minVal, maxVal, length=interVals)
  indexes=NULL
  sOrder=NULL
  maxi_densi=NULL
  densi_aux=densi
  
  # den_5000 <- c(0.00382540052, 0.05644854149, 0.31424556957, 0.84868620578,
  #              1.35600632563, 1.49691078926, 1.28723231246, 0.94488526689, 0.63404109120,
  #               0.40890113163, 0.26256638550, 0.17190027246, 0.11642706784, 0.08218805420,
  #               0.06059187920,0.04655503531, 0.03708990838, 0.03041913886, 0.02546148620,
  #               0.02154413757,0.01824544249, 0.01531040853, 0.01260439781, 0.01008346765,
  #               0.00776811870,0.00571433775, 0.00398268511, 0.00261137875, 0.00160061956,
  #               0.00091198894)
  
  den_5000 <- rep(0, 5000)
  
  #average of the reconstructions
  aver=NULL
  
  #max of each density
  maxi_densi=apply(densi, 1, max) 
  
  if (length(which(maxi_densi==0))!=0){
    densi=densi[-which(maxi_densi==0),]
    maxi_densi=maxi_densi[-which(maxi_densi==0)]
  }
  sOrder=sort.int(maxi_densi, index.return = TRUE)
  
  aver=apply(densi, 2, mean)  #without zeros
  aver2=apply(densi_aux, 2, mean) #with zeros
  
  #linea sup and inf (smaller first, greater last)
  #lin_inf=apply(densi[sOrder$ix[1:100],], 2, min) #densi[sOrder$ix[1],] 
  #lin_sup=apply(densi[sOrder$ix[(length(maxi_densi)-100):length(maxi_densi)],], 2, max) #densi[sOrder$ix[length(maxi_densi)],] 
  lin_inf=apply(densi, 2, min) #densi[sOrder$ix[1],] 
  lin_sup=apply(densi, 2, max) #densi[sOrder$ix[length(maxi_densi)],] 
  
  indexes=c(sOrder$ix[1:100],sOrder$ix[(length(maxi_densi)-100):length(maxi_densi)])
  
  x11()
  #histo=hist(S,freq = FALSE,breaks=5,main="",ylab="density",xlab="Severity",ylim=c(0,6)) #br=5
  plot(x,densi[1, ],type="l",pch=23,col="gray",lty=1,lwd=3,ylim=c(0,max(maxi_densi)+0.1),
       main="",ylab="density",xlab="Losses", xlim=c(minVal,maxVal+0.1))
  grid(col = "lightgray", lty = "dotted")
  
  polygon(c(x,rev(x)),c(densi[indexes[length(indexes)],],rev(densi[indexes[1],])),col="lightgray",border=NA)  
  
  for (i in 1:(length(indexes)/2)) 
  {
    polygon(c(x,rev(x)),c(densi[indexes[length(indexes)-i],],rev(densi[indexes[i],])),col="lightgray",border=NA)   
    #lines(spline(x,densi[indexes[i],]),type="l",pch=23,col="gray87",lty=1,lwd=3) #"gray87" "lightgray"
  }
  
  lines(spline(x, aver))
  #lines(spline(x, den_5000), col = 'red')
  legend((max(x)-min(x))*porc,max(maxi_densi)+0.15, c("True","AVERAGE","Reconstructions"), lty=c(1,1,1),
         lwd=c(3,1,3),col=c("black","gray30","gray87"),pch=c(NA,19,NA),bg = "white") 

  return(0)
}
#------------------------------------------------------------------------------------------------------------------------


# GRAPH
simGraph_x2 <- function(densi, S, minVal=0, maxVal=2, X2, porc=1.83, density_TRUE = 1){
  X2         <- as.vector(as.matrix(X2))
  interVals  <- dim(densi)[2]
  x          <- seq(minVal, maxVal, length = interVals)
  indexes    <- NULL
  sOrder     <- NULL
  maxi_densi <- NULL
  densi_aux  <- densi
  
  
  if (density_TRUE == 1){ # paint the TRUE density
    den_5000 <- c(0.00382540052, 0.05644854149, 0.31424556957, 0.84868620578,
                1.35600632563, 1.49691078926, 1.28723231246, 0.94488526689, 0.63404109120,
                0.40890113163, 0.26256638550, 0.17190027246, 0.11642706784, 0.08218805420,
                0.06059187920,0.04655503531, 0.03708990838, 0.03041913886, 0.02546148620,
                0.02154413757,0.01824544249, 0.01531040853, 0.01260439781, 0.01008346765,
                0.00776811870,0.00571433775, 0.00398268511, 0.00261137875, 0.00160061956,
                0.00091198894)
  }else{   # maxent density
    alpha       <- c(1.5000000, 0.7500000, 0.5000000, 0.3750000, 0.3000000,
                      0.2500000, 0.2142857, 0.1875000)
    lambda_5000 <- c(16240.12132,  -748.02989, -5997.50028,  1105.38763,  
                     4745.65939,  3727.93062, -149.98752, -5356.22172)
    x            <- sort(S[S>0])
    F_SME        <- F_diff_new(x, lambda_5000, alpha, globalZ = 1e-5)
    den_5000     <- F_SME$F_emp
  }
  
  #average of the reconstructions
  aver <- NULL
  #max of each density
  maxi_densi <- apply(densi, 1, max) 
  
  if (length(which(maxi_densi == 0)) != 0){
    densi      <- densi[-which(maxi_densi == 0),]
    maxi_densi <- maxi_densi[-which(maxi_densi == 0)]
  }
  sOrder <- sort.int(maxi_densi, index.return = TRUE)
  
  #averagw
  aver  <- apply(densi,     2, mean)  #without zeros
  aver2 <- apply(densi_aux, 2, mean) #with zeros
  
  #linea sup and inf (smaller first, greater last)
  #lin_inf=apply(densi[sOrder$ix[1:100],], 2, min) #densi[sOrder$ix[1],] 
  #lin_sup=apply(densi[sOrder$ix[(length(maxi_densi)-100):length(maxi_densi)],], 2, max) #densi[sOrder$ix[length(maxi_densi)],] 
  lin_inf <- apply(densi, 2, min) #densi[sOrder$ix[1],] 
  lin_sup <- apply(densi, 2, max) #densi[sOrder$ix[length(maxi_densi)],] 
  
  indexes <- c(sOrder$ix[1:100],sOrder$ix[(length(maxi_densi)-100):length(maxi_densi)])
  
  png(paste0("figures/myplot", TAM, ".png"), width=6, height=6, units="in", res=300)
  #x11()
  #histo=hist(S,freq = FALSE,breaks=5,main="",ylab="density",xlab="Severity",ylim=c(0,6)) #br=5
  plot(X2,densi[1, ],type="l",pch=23,col="gray",lty=1,lwd=3,ylim=c(0,max(maxi_densi)+0.1),
       main="",ylab="density",xlab="Losses", xlim=c(minVal,maxVal+0.1))
  grid(col = "lightgray", lty = "dotted")
  
  polygon(c(X2,rev(X2)), c(densi[indexes[length(indexes)],], rev(densi[indexes[1],])), col="lightgray", border=NA)  
  
  for (i in 1:(length(indexes)/2)) {
    polygon(c(X2,rev(X2)),c(densi[indexes[length(indexes)-i],],rev(densi[indexes[i],])),col="lightgray",border=NA)   
    #lines(spline(x,densi[indexes[i],]),type="l",pch=23,col="gray87",lty=1,lwd=3) #"gray87" "lightgray"
  }
  
  lines(spline(X2, aver), lwd=1, col = "gray30")
  if (density_TRUE == 1){
    x=seq(minVal, maxVal, length=length(den_5000))
    lines(spline(x, den_5000), col = 'red', lwd=3)
  }else {
    x=seq(minVal, maxVal, length=length(den_5000))
    lines(x, den_5000, col = 'red', lwd=3)
    #x11()
    print(class(den_5000))
    #plot(x, den_5000, col = 'red', lwd=3)
  }
  legend((max(x)-min(x))*porc,max(maxi_densi)+0.15, c("True", "AVERAGE", "Reconstructions"), lty=c(1,1,1),
#         lwd=c(3,1,3),col=c("red","gray30","gray87"), pch=c(NA,19,NA),bg = "white") 
          lwd=c(3,1,3),col=c("red","gray30","gray87"),bg = "white") 
 dev.off()
  return(0)
}

