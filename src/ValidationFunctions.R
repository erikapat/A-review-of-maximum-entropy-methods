
# metrics to validate the results
metrics <- function(F_maxent, F_e){
  
  #MEAN ERROR
  ME <- (1/length(F_maxent))*sum(F_maxent-F_e)
  #MEAN ABSOLUTE ERROR #L1-norm
  MAE <- (1/length(F_maxent))*sum(abs(F_maxent-F_e))
  #MEAN SQUARED ERROR (#BRIER SCORE: measures the mean squared error. Negative orientation (smaller score the better))
  MSE <- (1/length(F_maxent))*sum((F_maxent-F_e)^2)
  #L2-norm
  #ROOT MEAN SQUARED ERROR
  RMSE <- sqrt((1/length(F_maxent))*sum((F_maxent-F_e)^2))
  
  metrics <- data.frame(ME = ME, MAE = MAE, MSE = MSE, RMSE = RMSE)
  
  return(metrics)
}


PIT <- function(F_maxent, perc = .1){
  #PIT
  
  if (length(F_maxent) < 50){
    perc = 1
  }
  sampF <- sample(F_maxent, length(F_maxent)*perc, replace = FALSE, prob = NULL)
  sampF <- data.frame(sampF = sampF)
  gg <- ggplot(sampF, aes(x = sampF)) + geom_histogram(bins = sqrt(length(sampF$sampF)) )
  gg <- gg +  theme_tq() 
  gg <- gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10), axis.title.x = element_text(color="black", size= 12, face="bold"))
  gg <- gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12, face="bold"))
  gg <- gg + theme(legend.title = element_text(colour="blue", size=10, 
                                                           face="bold")) 
  gg <- gg +   guides(col = guide_legend(ncol = 1)) 
  gg <- gg + ggtitle("PIT")   + xlab("CDF MAXENT") + ylab("density")  + theme(plot.title = element_text(size = 15, face = "bold"))
  
 return(gg)
  
}

#------------------------

tests <- function(F_maxent){
  #test
  
  ##############################
  #    TESTS
  ##############################
  
  if (size_test >2000) #I assume that we use validation data
  {
    cat("\n","\t","ERROR: size_test should be less than 2000","\n")
    return(0)
  }
  if (size_test <=0) #I assume that we use validation data
  {
    cat("\n","\t","ERROR: size_test should greater than 0","\n")
    return(0)
  }
  
  aux_F=F_maxent[F_maxent>0]
  aux_F=aux_F[aux_F<1]
  F_aux1=sample(aux_F, length(aux_F), replace = FALSE, prob = NULL)
  #if (length(aux_F)>2000)
  #{
  #  F_aux1=NULL
  #  F_aux1=sample(aux_F, size_test, replace = FALSE, prob = NULL)
  #}
  
  NN=length(F_aux1) #F_aux1
  cat(rep("\n", 1)) 
  cat("TESTS:","\n")
  cat("-----------------------------------------------------------","\n")
  
  #KOLMOGOROV TEST
  cat("KOLMOGOROV TEST:","\n")
  pksR<-ks.test(F_aux1, punif,alternative = "two.sided")
  cat("p.value: ",pksR$p.value,"\n") 
  cat("statistic: ",sqrt(NN)*pksR$statistic,"\n") 
  cat("Critical values: ","\n")
  cat(rep("\t",2),"1%: ",1.63,"\n")
  cat(rep("\t",2),"5%: ",1.36,"\n")
  #cat(rep("\t",2),"1%: ",sqrt(length(F))*unif.test.quantile(type= "ks", length(F), 0.01),"\n")
  #cat(rep("\t",2),"5%: ",sqrt(length(F))*unif.test.quantile(type= "ks", length(F), 0.05),"\n")
  #cat("KOLMOGOROV TEST (simulated):","\n")
  #pKS=p_value_KS(F,rep=1000,sample_correction=0) #sqrt(N)*D_n
  #cat("p.value: ",pKS$p_value,"\n") 
  #cat("statistic: ",pKS$statistic,"\n") 
  #cat("Critical values: ","\n")
  #cat(rep("\t",2),"1%: ",pKS$critical_value_001,"\n")
  #cat(rep("\t",2),"5%: ",pKS$critical_value_005,"\n")
  
  #aux_F=F[F>0]
  #aux_F=aux_F[aux_F<1]
  cat(rep("\n", 1)) 
  cat("-----------------------------------------------------------","\n")
  cat("ANDERSON-DARLING TEST: ","\n") 
  cat("p.value: ",ad.test(F_aux1, punif)$p.value,"\n")
  cat("Statistic: ",ad.test(F_aux1, punif)$statistic,"\n")
  #cat("ANDERSON-DARLING TEST: (simulated):","\n")
  #An=p_value_anderson(aux_F,rep=3000)
  #cat("p.value: ",An$p_value,"\n")
  #cat("Statistic: ",An$statistic,"\n")
  #cat("Critical values: ","\n")
  #cat(rep("\t",2),"1%: ",An$critical_value001,"\n")
  #cat(rep("\t",2),"5%: ",An$critical_value005,"\n")
  cat("Critical values: ","\n")
  cat(rep("\t",2),"1%: ",3.857,"\n")
  cat(rep("\t",2),"5%: ",2.492,"\n")
  
  
  cat(rep("\n", 1)) 
  cat("-----------------------------------------------------------","\n")
  cat("CRAMER-VON MISES test: ","\n")
  cat("Statistic: ", unif.test.statistic(F_aux1,type= "cvm"),"\n")
  #critical values
  cat("Critical values: ","\n")
  cat(rep("\t",2),"1%: ",unif.test.quantile(type= "cvm", length(F), 0.01),"\n")
  cat(rep("\t",2),"5%: ",unif.test.quantile(type= "cvm", length(F), 0.05),"\n")
  #cat("CRAMER-VON MISES test (simulated): ","\n")
  #Wn=p_value_cramer(F,rep=1000)
  #cat("p.value: ",Wn$p_value,"\n")
  #cat("Statistic: ",Wn$statistic,"\n")
  #cat("Critical values: ","\n")
  #cat(rep("\t",2),"1%: ",Wn$critical_value001,"\n")
  #cat(rep("\t",2),"5%: ",Wn$critical_value005,"\n")
  
  cat(rep("\n", 1)) 
  cat("-----------------------------------------------------------","\n")
  
  #Transform data to a Normal N(0,1).
  #sampF_aux=F[F<1] #*** 1's generates Inf
  #sampF_aux=sampF_aux[sampF_aux>0] #*** 0's generates -Inf
  nvector = qnorm(F_aux1)
  
  if (length(F[F==1])>0) 
  {
    cat("Warning: INF's in Normal transformation (qnorm(1)=Inf)","\n")
    #output
    #return(list(ME=ME,MAE=MAE,MSE=MSE,RMSE=RMSE,log_score=log_score))
  }
  
  nvector_aux=sample(nvector,length(nvector), replace = FALSE, prob = NULL)
  
  #x11()
  #hist(nvector_aux,main="Normal Transformation",xlab="Normal Values")
  
  
  cat(rep("\n", 1)) 
  cat("-----------------------------------------------------------","\n")
  cat("Berkowitz Test:","\n")
  cat("\t","(1) LR1 -> Normal (0,1):","\n")
  Ber=BerkowitzTest(data = nvector, lags = 1, significance = 0.01)
  cat(rep("\t", 2),"pvalue: ",Ber$LRp,"\n")
  cat(rep("\t", 2),"statistic: ",Ber$LR,"\n")
  cat(rep("\t", 2),"NULL hypothesis: ",Ber$H0,"\n")
  cat(rep("\t", 2),"Decision: ",Ber$Decision," at significance = 0.01","\n")
  cat("Critical values: ","\n")
  cat(rep("\t",2),"1%: ",qchisq(0.99,3),"\n")
  cat(rep("\t",2),"5%: ",qchisq(0.95,3),"\n")
  
  #this is for the tail
  #cat("\t","(2) LR3 -> Normal (0,1) with autocorrelation:","\n") 
  #Normal with autocorrelation LR1
  #Ber2=BerkowitzTest(data = nvector, significance = 0.01, tail.test=TRUE)
  #cat(rep("\t", 2),"pvalue: ",Ber2$LRp,"\n")
  #cat(rep("\t", 2),"statistic: ",Ber2$LR,"\n")
  #cat(rep("\t", 2),"NULL hypothesis: ",Ber2$H0,"\n")
  #cat(rep("\t", 2),"Decision: ",Ber2$Decision," at significance = 0.01","\n")
  
  #cat("Critical values: ","\n")
  #cat(rep("\t",2),"1%: ",qchisq(0.99,3),"\n")
  #cat(rep("\t",2),"5%: ",qchisq(0.95,3),"\n")
  #cat("-----------------------------------------------------------","\n")
  #cat(rep("\n", 1))
  cat("-----------------------------------------------------------","\n")
  
  #set.seed(1)
  #if (length(nvector)>4000)
  #{
  #  nvector=sample(nvector, length(nvector), replace = FALSE, prob = NULL)
  #}
  #cat("JARQUE-BERA TEST: ","\n")
  #cat("p.value: ",jarque.bera.test(nvector_aux)$p.value,"\n")
  #cat("statistic: ",jarque.bera.test(nvector_aux)$statistic,"\n")
  cat("ROBUST JARQUE-BERA TEST : ","\n")
  stJB=rjb.test(nvector_aux,option = "RJB",crit.values = "chisq.approximation")
  cat("p.value: ",stJB$p.value,"\n")
  cat("statistic: ",stJB$statistic,"\n")
  stJB=rjb.test(nvector,option = "RJB",crit.values = "empirical", N = 1000)
  cat("p.value: ",stJB$p.value,"\n")
  cat("statistic: ",stJB$statistic,"\n")
  cat("JARQUE-BERA TEST (simulated): ","\n")
  JPn=p_value_JP(nvector,rep=1000)
  cat("p.value: ",JPn$p_value,"\n")
  cat("statistic: ",JPn$statistic,"\n")
  cat("Critical values (Asymptotic): ","\n")
  #cat(rep("\t",2),"1%: ",JPn$critical_value001,"\n")
  #cat(rep("\t",2),"5%: ",JPn$critical_value005,"\n")
  cat(rep("\t",2),"1%: ",qchisq(0.99,2),"\n")
  cat(rep("\t",2),"5%: ",qchisq(0.95,2),"\n")
  
  
  ##############
  #independence#
  ##############
  sampF=sample(F, length(F), replace = FALSE, prob = NULL)
  sampF_e=sample(F_e, length(F_e), replace = FALSE, prob = NULL)
  #sampF=sample(F, size_test, replace = FALSE, prob = NULL)
  #sampF_e=sample(F_e, size_test, replace = FALSE, prob = NULL)
  
  x11()
  par(mfrow=c(2,2)) 
  ll=hist(sampF, breaks=5, main="PIT histogram")
  acf(sampF, main="ACF of PIT") #ACF OF PIT
  acf(sampF-1/2, main="ACF of PIT-1/2")
  acf((sampF-1/2)^2, main="ACF of (PIT-1/2)^2")
  #acf((sampF-1/2)^3)
  #partial
  x11()
  par(mfrow=c(2,2)) 
  pacf(sampF, main="PACF of PIT")
  pacf(sampF-1/2, main="ACF of PIT-1/2")
  pacf((sampF-1/2)^2, main="ACF of (PIT-1/2)^2")
  pacf((sampF-1/2)^3, main="ACF of (PIT-1/2)^3")
  #######################################################################
}

