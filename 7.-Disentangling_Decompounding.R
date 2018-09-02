

########################################################################
#SIMULATIONS

rm(list=ls()) #remover cualquier variable del espacio de trabajo :)
require(ggplot2)
require(tidyquant)
require(data.table)
source("src/disentangling.R")
dir <- paste0('data/', 'frequency_two_source_risk.dat')
n = fread(dir)

#-------------------------------------------------------------------------------------------------------------------------------------

freq <- as.data.table(n)
names(freq) <- c('n')
densi_gg <- ggplot(freq) + geom_histogram(aes(x = n, y = ..density..), bins = 10, color = 'gray', fill = 'darkgray')  # 
densi_gg <- densi_gg +  theme_tq() 
densi_gg <- densi_gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12, face="bold"))
densi_gg <- densi_gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12, face="bold"))
densi_gg <- densi_gg + theme(legend.title = element_text(colour="blue", size=10, 
                                                         face="bold")) 
densi_gg <- densi_gg +   guides(col = guide_legend(ncol = 1)) 
densi_gg <- densi_gg + ggtitle("FREQUENCY DISTRIBUTION")   + xlab("Frequency of Losses") + ylab("density")  + theme(plot.title = element_text(size = 15, face = "bold"))

dir <- '/home/erika/Dropbox/4ErikaIII (continuacion)/(R) CODIGO R MAXENT NEW/1.-TEXTO_PAPER/NUMERICAL_CASE_PAPER/'
if(!is.null(dir)) ggsave(paste0(dir, '/FREQ_DIST', '.png'))
#-------------------------------------------------------------------------------------------------------------------------------

nn <- dis_data_improve(n$V1) 
nn <- nn[2:nrow(nn), ]
###############################
matr1=recursion(nn)
matr1=round(matr1)

list=EM_method(matr1,mod=2,G=2)
xx <- do_groups(matr1,list$list,dist="Poisson",lineas=0,kmeans=1)

#--------------------------------------------------------
# GRAPH 
clusters=kmeans(matr1[,-2],2,iter.max=100,nstart=25)
points(clusters$centers, pch = c(8,4), cex=5)

dat_freq <- as.data.table(cbind(matr1, clusters$cluster))
names(dat_freq) <- c('index', 'freq', 'akb', 'cluster')
dat_freq <- dat_freq[-1, ]

centers <- as.data.frame(clusters$centers)
gg <- ggplot(dat_freq) + geom_point(aes(x = index, y = akb, color = factor(cluster)), size = 5,fill = 'black',  show.legend=F)  # 
gg <- gg +  theme_tq() + ylim(c(-10,20))
gg <- gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12, face="bold"))
gg <- gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12, face="bold"))
gg <- gg + theme(legend.title = element_text(colour="blue", size=10, 
                                                         face="bold")) 
gg <- gg + scale_color_manual(name="", values = c("black", "gray"))
#gg <- gg +   guides(col = guide_legend(ncol = 1)) 
gg <- gg + ggtitle("Panjer Recursion Formula")   + xlab("Frequency of Losses") + ylab("recursion values")  + theme(plot.title = element_text(size = 15, face = "bold"))
gg

gg <- gg + geom_point(data=centers, aes(x=V1,y=V2), shape = 1, size = 8) 
gg <- gg + geom_point(data=centers, aes(x=V1,y=V2), shape = 5, size = 8) 
gg <- gg + geom_point(data=centers, aes(x=V1,y=V2), shape = 3, size = 8) 
#gg <- gg + geom_point(data=centers, aes(x=V1,y=V2), size=52, alpha=.3, legend=FALSE)
gg
dir <- '/home/erika/Dropbox/4ErikaIII (continuacion)/(R) CODIGO R MAXENT NEW/1.-TEXTO_PAPER/NUMERICAL_CASE_PAPER/'
if(!is.null(dir)) ggsave(paste0(dir, '/PanjerRecursion', '.png'))


#-------------------------------------------------------------------------------------------------------------------------------------------------------

# #media por max verosimilitud
# lam_t=mean(c(v_group1,v_group2))
 lam=4 #mean(v_group1)  
 lam2=8 #mean(v_group2) 
# 
# freq_obs=matr1[,2] #frecuencia observada
# 
# 
# i=matr1[,1]
# p1=(lam_t-lam2)/(lam-lam2)
#pk_aux=p1*((exp(-lam)*lam^i)/factorial(i))+ (1-p1)*((exp(-lam2)*lam2^i)/factorial(i))
NN = 4 + 8
NN1 = 4 
NN2 =  8
i = 0:20
pk_aux=(NN1/NN)*((exp(-lam)*lam^i)/factorial(i))+ (NN2/NN)*((exp(-lam2)*lam2^i)/factorial(i))
dt = data.frame(x = i, densi = pk_aux)

freq <- as.data.table(n)
names(freq) <- c('n')
densi_gg <- ggplot(freq) + geom_histogram(aes(x = n, y = ..density..), bins = 10, color = 'gray', fill = 'darkgray')  # 
densi_gg <- densi_gg + geom_line(data = dt, aes(x = x, y = densi)) 
densi_gg <- densi_gg +  theme_tq() 
densi_gg <- densi_gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12, face="bold"))
densi_gg <- densi_gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12, face="bold"))
densi_gg <- densi_gg + theme(legend.title = element_text(colour="blue", size=10, 
                                                         face="bold")) 
densi_gg <- densi_gg +   guides(col = guide_legend(ncol = 1)) 
densi_gg <- densi_gg + ggtitle("FREQUENCY DISTRIBUTION")   + xlab("Frequency of Losses") + ylab("density")  + theme(plot.title = element_text(size = 15, face = "bold"))

dir <- '/home/erika/Dropbox/4ErikaIII (continuacion)/(R) CODIGO R MAXENT NEW/1.-TEXTO_PAPER/NUMERICAL_CASE_PAPER/'
if(!is.null(dir)) ggsave(paste0(dir, '/FREQ_DIST', '.png'))
#---------------------------------------------------------
# 
# source("disentangling.r")
# bin_parameters(matr1[,-2]) #regresion de toda la recta
# 
# x11()
# hist(c(n1,n2),breaks=28)
# barplot(n[,2],xlab = 'Number of errors', ylab = 'frequency', col='grey99')
# 
# n1=read.table(file ="n1_1.dat", header = FALSE)
# n2=read.table(file ="n2_1.dat", header = FALSE)
# histo1=hist(n1[,1])
# histo2=hist(n2[,1])
# x11()
# plot(histo2$mids,histo2$counts,type='s', 
#      ylim=c(0, max(histo1$counts, histo2$counts)),
#      xlim=c(0, max(histo1$mids, histo2$mids)),
#      xlab="n",ylab="frequencies")
# lines(histo1$mids,histo1$counts,type='s')
# 
# 
# lines(histo2$mids,histo2$counts,type='s')
# 
# ################################################
# 
# sh=5
# rt=15
# 
# dat_g=NULL
# S1=vector(mode = "numeric", length = M1)
# indi=NULL
# # for each k 
# for (k in 1:M1)
# {
#   xc=NULL
#   xc=(rgamma(n1[k], shape = sh, rate = rt)) 
#   S1[k]=sum(xc)
#   dat_g=c(dat_g,xc)
# }
# 
# #GRAPH of the simulation
# hist(S1[S1>0],breaks=15,freq=FALSE,include.lowest = TRUE)
# 
# sh=5
# rt=15
# 
# S2=vector(mode = "numeric", length = M2)
# indi=NULL
# # for each k 
# dat_g2=NULL
# for (k in 1:M2)
# {
#   xc=NULL
#   xc=(rgamma(n2[k], shape = sh, rate = rt)) 
#   S2[k]=sum(xc)
#   dat_g2=c(dat_g2,xc)
# }
# #GRAPH of the simulation
# hist(S2[S2>0],breaks=15,freq=FALSE,include.lowest = TRUE)
# 
# hist(c(dat_g2,dat_g),breaks=15,freq=FALSE)
# 
# S=NULL
# S2[(length(S2)+1):length(S1)]=0
# #S2[1:(length(S1)-(length(S2)+1))]=0
# #S=c(S2+S1[1:length(S2)],S1[(length(S2)+1):length(S1)])
# S=S1+S2
# hist(S[S>0],breaks=15,freq=FALSE,include.lowest = TRUE)
# 
# #save
# write.table(n1, file = "n1_1.dat",row.names = FALSE,col.names = FALSE)
# write.table(n2, file = "n2_1.dat",row.names = FALSE,col.names = FALSE)
# write.table(n, file = "frequencies1.dat",row.names = FALSE,col.names = FALSE)
# write.table(S, file = "prueba1.dat",row.names = FALSE,col.names = FALSE)
# write.table(matr1[,-2], file = "recursion1.dat",row.names = FALSE,col.names = FALSE)
# 
# 
# #source("general_functions.R")
# #simul_compound_gamma("prueba2.dat", l=8, a=2, b=3, 500)
# 
# #ALL above can be simplified with the function 
# #simul_compound_gamma<-function(name, ele, a=2, b=1, M)
# #in the file general_functions.R in the folder R code in ERIKAIII 
# ##############################################################################
# 
# 
# ##############
# #EMPEZAR AQUI#
# ##############
# 
# rm(list=ls()) #remover cualquier variable del espacio de trabajo :)
# setwd("C:\\Users\\Erika\\Dropbox\\4ErikaIII (continuacion)\\paper (4)\\code")
# 
# source("opt_methods.r")
# source("SME_functions.r")
# 
# set.seed(100)
# #cargar data
# SS=read.table(file ="prueba1.dat", header = FALSE)
# #transformar a vector
# S=SS[,1]
# dim(SS)
# k=8
# length(S)
# hist(S)
# 
# mp=maxent_parameters(S,k) #aqui se deben incluir los ceros (OJO)
# mu=mp$mu
# alpha=mp$alpha
# 
# 
# global<<-0.0001 #0.0055 #0.0001 #number of points in the integration (global variable)
# #shorting this value the convergence can be more fast.
# 
# #fix always a tolerance of 1e-2 gives good results in general, in some cases
# #you can improved, but not always. Depends of the case.
# tol=4.4e-7
# lambda=rep(1,k)
# source("opt_methods.r")
# #Initial values
# h=BB_modified(lambda,alpha,mu,M=5000,tolerance=tol,step=2,correction=0,GLL=1,
#               sigma1=0.1,sigma2=0.5,phi_max=1e20,phi_min=1e-20,N=10)
# if (h$error==000 | h$error==-999 | h$error==-666) error(h$error)
# lambda=h$lambda_Opt
# 
# #################
# # density
# #################
# #density and graphs of S and estimated
# source("SME_functions.r")
# x=seq(0,18,length=30) #x o S.
# den=SME_density(lambda,alpha,x,S,globalZ=1e-5)
# Z=den$Z
# fS=den$densidadS
# 
# source("SME_functions.R")
# ll=verification(S,lambda,alpha,gamma=0.92,size_test=1000,
#                 globalZ=1e-5) #0.0055
# ll
# 
# ########################
# #DECOMPOUNDING
# ########################
# 
# source("SME_functions.R")
# poisson_par=19
# hu=laplace_individuals(lambda,alpha,poisson_par,globalZ=1e-1,mu) #0.00155
# hu
# 
# #BB method...
# tol=1e-4
# lambda_p=rep(1,8)
# source("opt_methods.r")
# hh=BB_modified(lambda_p,alpha,hu,M=1000,tolerance=tol,step=2,correction=0,GLL=1,
#                sigma1=0.1,sigma2=0.5,phi_max=1e20,phi_min=1e-20,N=10)
# if (h$error==000 | h$error==-999 | h$error==-666) error(h$error)
# lambda_p=hh$lambda_Opt
# 
# #density and graphs of S and estimated
# source("SME_functions.r")
# x=seq(0,1,length=50) #x o S.
# den=SME_density(lambda_p,alpha,x,S=NULL,globalZ=1e-5)
# Z=den$Z
# fS=den$densidadS
# plot(x,fS,type="o")
# 
# 
# x11()
# #cargar data
# source("disentangling.r")
# indi=mix_simul_gamma(size1=400,size2=500, seed=1)
# hist(indi,freq = FALSE,breaks=15,main="",ylab="density",xlab="Severity",)
# simullog=4/19*dgamma(x,shape=5,rate=15) +16/19*dgamma(x,shape=5,rate=15)
# lines(x,simullog,ylim = c(0,max(simullog,fS)+0.1),main="",ylab="density",xlab="Severity",col="black",type="l",lwd=3.5)
# lines(x,fS,type="l",col="gray",lty=2,lwd=3.5)
# legend(max(x)-2,2, c("Real","SME"), lty=c(1,4), lwd=c(3.5,2.5),col=c("black","gray")) 
# 
# #OTRO CASO
# #########################################################################
# ########################################################################
# #SIMULATIONS (EJEMPLO PAPER !!!!!!!)
# 
# rm(list=ls()) #remover cualquier variable del espacio de trabajo :)
# setwd("C:\\Users\\Erika\\Dropbox\\4ErikaIII (continuacion)\\paper (4)\\code")
# 
# #simulation of Poiss-Gamma
# #1.-Equal sizes, all is independent, different poisson parameters, 
# #same gamma parameters
# 
# #sizes
# M1=500
# M2=500
# #parameter
# ele1=4
# ele2=6
# 
# #M data from a poisson of lambda =1, k=1,...M
# 
# #dist Poisson
# #seed=5
# #set.seed(seed)
# n1=rpois(M1,ele1) 
# n2=rpois(M2,ele2)
# 
# source("disentangling.r")
# n=dis_data(n1,n2)
# 
# ###############################
# source("disentangling.r")
# n=read.table(file ="frequencies6.dat", header = FALSE)
# n_r=n
# matr1=recursion(n)
# #matr1=round(matr1)
# #matr1=matr1[c(-14),]
# list=EM_method(matr1,mod=1,G=3)
# do_groups(matr1,list$list,dist="Poisson",lineas=0,kmeans=1,ylim2=8)
# #points(matr1[14,1],matr1[14,3],pch=19)
# 
# 
# #CRITERIOS PARA ESTIMAR EL NUMERO DE PARAMETROS
# #1: By rule of thumb
# sqrt((dim(matr1)[1])/2)
# #2: Elbow method
# mydata<-matr1[,-2]
# # Determine number of clusters
# wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# for (i in 2:13) wss[i] <- sum(kmeans(mydata, 
#                                      centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares")
# 
# #3: Information Criterion Approach
# source("disentangling.r")
# Grupos=4
# list1=EM_method(matr1,mod=3,G=Grupos)
# d=list1$degrees_of_freedom
# -2*list1$loglikelihood
# list1$bic
# list1$BIC
# num_grupos(log_lik=list1$loglikelihood,ENT=list1$entropy, grupos=Grupos,n=dim(n_r)[1],printi=1,degrees_of_freedom=d) 
# list1$negentropy  
# 
# #################################################################################
# 
# ################################################################################
# #### PRUEBAS DE HIPOTESIS
# 
# lis=list$list
# group1=NULL
# group2=NULL
# group1=matr1[lis==1,]
# group2=matr1[lis==2,]
# 
# v_group1=rep(group1[,1],group1[,2])
# sum(group1[,2])
# length(v_group1)
# mean(v_group1)
# var(v_group1)
# skewness(v_group1) 
# 
# v_group2=rep(group2[,1],group2[,2])
# mean(v_group2)
# var(v_group2)
# skewness(v_group2,type=3) 
# 
# NN1=sum(matr1[lis==1,2])
# NN2=sum(matr1[lis==2,2])
# NN=NN1+NN2
# 
# NN1/NN
# NN2/NN
# 
# 
# #################################  
# #media por max verosimilitud
# lam_t=mean(c(v_group1,v_group2))
# lam=4 #mean(v_group1)  
# lam2=6 #mean(v_group2) 
# 
# freq_obs=matr1[,2] #frecuencia observada
# 
# 
# i=matr1[,1]
# p1=(lam_t-lam2)/(lam-lam2)
# #pk_aux=p1*((exp(-lam)*lam^i)/factorial(i))+ (1-p1)*((exp(-lam2)*lam2^i)/factorial(i))
# pk_aux=(NN1/NN)*((exp(-lam)*lam^i)/factorial(i))+ (NN2/NN)*((exp(-lam2)*lam2^i)/factorial(i))
# 
# sum(pk_aux)
# e=NN*pk_aux #frecuencia teorica
# cbind(freq_obs,e) #tabla, combinar los menores que 5
# freq_obs1=c(freq_obs[1:11],sum(freq_obs[e<=8]))
# e1=c(e[1:11],sum(e[e<=8]))
# cbind(freq_obs1,e1)
# #e1=e
# #freq_obs1=freq_obs
# 
# val_chi=((freq_obs1-e1)^2)/(e1)
# chi=sum(  val_chi  )
# chi
# df=length(val_chi)-1
# qchisq(0.999, df)
# #no reject H0
# #p-value
# 1-pchisq(chi,df) #if p_value if grater than alpha, we dont reject Ho
# 
# 
# ######################GRAFICO ###################################################################3
# i=0:30
# p1=(lam_t-lam2)/(lam-lam2)
# #pk_aux=p1*((exp(-lam)*lam^i)/factorial(i))+ (1-p1)*((exp(-lam2)*lam2^i)/factorial(i))
# pk_aux=(NN1/NN)*((exp(-lam)*lam^i)/factorial(i))+ (NN2/NN)*((exp(-lam2)*lam2^i)/factorial(i))
# ke=NN*pk_aux
# 
# 
# x11()
# barplot(n[,2],xlab = 'Number of Losses', ylab = 'frequency', col='grey99',ylim=c(0,190))
# xxx=c(1:length(ke))
# lines(cbind(c(0,xxx+1),c(0,ke)),lwd=3) 
# lines(cbind(c( 14.5,15.5,16.5,17.5,18.5,19,20.5),c(0,0,0,0,0,0,0)),lwd=4)
# 
# 
# 
# #############################################################################
# #Poisson plot 
# pos_plot=log(n[-c(14,15),2])+log(factorial(n[-c(14,15),1]))
# plot(n[-c(14,15),1],pos_plot)
# 
# mc=cbind(n[-c(14,15),1],n[-c(14,15),1],pos_plot)
# 
# list1=EM_method(mc,mod=4,G=2)
# list=list1$lista
# hj=do_groups(mc,list,dist="Poisson",lineas=1,kmeans=0,ylim2=40,
#              true0=1,xlabeli='k',ylabeli='log(xk)+log(k!)')
# hj$b_coef
# exp(hj$b_coef)
# 
# ###################################################################################################
# 
# 
# list1=EM_method(matr1,mod=3,G=1)
# num_grupos(log_lik=list1$loglikelihood,ENT=0, grupos=1,n=dim(n_r)[1],printi=1) 
# list1=EM_method(matr1,mod=3,G=2)
# num_grupos(log_lik=list1$loglikelihood,ENT=0, grupos=2,n=dim(n_r)[1],printi=1)
# list1=EM_method(matr1,mod=3,G=3)
# num_grupos(log_lik=list1$loglikelihood,ENT=0, grupos=3,n=dim(n_r)[1],printi=1)
# list1=EM_method(matr1,mod=3,G=4)
# num_grupos(log_lik=list1$loglikelihood,ENT=0, grupos=4,n=dim(n_r)[1],printi=1)
# list1=EM_method(matr1,mod=3,G=5)
# num_grupos(log_lik=list1$loglikelihood,ENT=0, grupos=5,n=dim(n_r)[1],printi=1)
# list1=EM_method(matr1,mod=3,G=6)
# num_grupos(log_lik=list1$loglikelihood,ENT=0, grupos=6,n=dim(n_r)[1],printi=1)
# 
# 
# 
# 
# 
# source("disentangling.r")
# bin_parameters(matr1[,-2]) #regresion de toda la recta
# 
# n1=read.table(file ="n1_6.dat", header = FALSE)
# n2=read.table(file ="n2_6.dat", header = FALSE)
# n1=n1[,1]
# n2=n2[,1]
# 
# x11()
# hist(c(n1[,1],n2[,1]),breaks=28)
# barplot(n[,2],xlab = 'Number of errors', ylab = 'frequency', col='grey99')
# 
# histo1=hist(n1[,1])
# histo2=hist(n2[,1])
# x11()
# plot(histo2$mids,histo2$counts,type='s', 
#      ylim=c(0, max(histo1$counts, histo2$counts)),
#      xlim=c(0, max(histo1$mids, histo2$mids)),
#      xlab="n",ylab="frequencies")
# lines(histo1$mids,histo1$counts,type='s')
# 
# 
# lines(histo2$mids,histo2$counts,type='s')
# 
# ################################################
# 
# sh=5
# rt=15
# 
# dat_g=NULL
# S1=vector(mode = "numeric", length = M1)
# indi=NULL
# # for each k 
# for (k in 1:M1)
# {
#   xc=NULL
#   xc=rgamma(n1[k], shape = sh, rate = rt)
#   S1[k]=sum(xc)
#   dat_g=c(dat_g,xc)
# }
# 
# #GRAPH of the simulation
# hist(S1[S1>0],breaks=15,freq=FALSE,include.lowest = TRUE)
# 
# #sh=2
# #rt=10
# 
# S2=vector(mode = "numeric", length = M2)
# indi=NULL
# # for each k 
# dat_g2=NULL
# for (k in 1:M2)
# {
#   xc=NULL
#   xc=(rgamma(n2[k], shape = sh, rate = rt)) 
#   S2[k]=sum(xc)
#   dat_g2=c(dat_g2,xc)
# }
# #GRAPH of the simulation
# hist(S2[S2>0],breaks=15,freq=FALSE,include.lowest = TRUE)
# x11()
# hist(c(dat_g2,dat_g),breaks=15,freq=FALSE)
# 
# S=NULL
# #S2[(length(S2)+1):length(S1)]=0
# #S2[1:(length(S1)-(length(S2)+1))]=0
# #S=c(S2+S1[1:length(S2)],S1[(length(S2)+1):length(S1)])
# S=S1+S2
# hist(S[S>0],breaks=15,freq=FALSE,include.lowest = TRUE)
# 
# #save
# write.table(n1, file = "n1_6.dat",row.names = FALSE,col.names = FALSE)
# write.table(n2, file = "n2_6.dat",row.names = FALSE,col.names = FALSE)
# write.table(n, file = "frequencies6.dat",row.names = FALSE,col.names = FALSE)
# write.table(S, file = "prueba6.dat",row.names = FALSE,col.names = FALSE)
# write.table(matr1[,-2], file = "recursion6.dat",row.names = FALSE,col.names = FALSE)
# write.table(c(dat_g2,dat_g), file = "individuals6.dat",row.names = FALSE,col.names = FALSE)
# 
# 
# 
# #source("general_functions.R")
# #simul_compound_gamma("prueba2.dat", l=8, a=2, b=3, 500)
# 
# #ALL above can be simplified with the function 
# #simul_compound_gamma<-function(name, ele, a=2, b=1, M)
# #in the file general_functions.R in the folder R code in ERIKAIII 
# ##############################################################################
# 
# 
# ##############
# #EMPEZAR AQUI#
# ##############
# 
# rm(list=ls()) #remover cualquier variable del espacio de trabajo :)
# setwd("C:\\Users\\Erika\\Dropbox\\4ErikaIII (continuacion)\\paper (4)\\code")
# 
# source("opt_methods.r")
# source("SME_functions.r")
# 
# set.seed(100)
# #cargar data
# SS=read.table(file ="prueba6.dat", header = FALSE)
# #transformar a vector
# S=SS[,1]
# dim(SS)
# k=8
# length(S)
# hist(S)
# 
# mp=maxent_parameters(S,k) #aqui se deben incluir los ceros (OJO)
# mu=mp$mu
# alpha=mp$alpha
# 
# 
# global<<-0.0001 #0.0055 #0.0001 #number of points in the integration (global variable)
# #shorting this value the convergence can be more fast.
# 
# #fix always a tolerance of 1e-2 gives good results in general, in some cases
# #you can improved, but not always. Depends of the case.
# tol=4.44e-7 # 4.44e-7 da el menor error (verificado ya)
# lambda=rep(1,8)
# source("opt_methods.r")
# #Initial values
# h=BB_modified(lambda,alpha,mu,M=10000,tolerance=tol,step=2,correction=0,GLL=1,
#               sigma1=0.1,sigma2=0.5,phi_max=1e20,phi_min=1e-20,N=10)
# if (h$error==000 | h$error==-999 | h$error==-666) error(h$error)
# lambda=h$lambda_Opt
# 
# #lambda=c(  36.9769199 -186.0850064  317.3737244  -30.3222614 -128.9302233 -77.4080173    6.5332488   72.5725473)
# #################
# # density
# #################
# graphics.off()
# #density and graphs of S and estimated
# source("SME_functions.r")
# x=seq(0,15,length=50) #x o S.
# den=SME_density(lambda,alpha,x,S,globalZ=1e-5)
# Z=den$Z
# fS=den$densidadS
# 
# densidad=fS
# 
# 
# source("SME_functions.R")
# ll=verification(S,lambda,alpha,gamma=0.92,size_test=1000,
#                 globalZ=1e-5) #0.0055
# ll
# ##################################ANALITYCAL DISTRIBUTION #############################################
# #sub_gamma(10,5,15,x,20)
# 
# #Poisson-gamma: Standard approach
# source("standard_approach.R")
# d12=sa_gamma(10, 5, 15,x,trap=1,globalZ=1e-3,maxi_n=25)
# 
# d12
# #nn=read.table(file ="frequencies6.dat", header = FALSE)
# #maxi_n_calculation(10,nn[,2],eps=1e-1)
# 
# ############################################################################
# # adaptar grafico GRAFICO
# sep1=15
# d=hist(S[S>0],breaks=sep1,freq=FALSE,include.lowest = TRUE,ylim = c(0,0.6),
#        xlab="S", main = paste("")) #main = paste("20 bars")
# 
# #PAPER
# x11()
# sep=15 #set 9 o 10
# par(ps = 12, cex.lab = 1.5, cex.main = 1.5) #, cex.main = .5,cex=.5,cex.lab=.5, cex.axis=.5, cex.main=.5
# hist(S[S>0],breaks=sep,freq=FALSE,include.lowest = TRUE,
#      ylim = c(0,max(c(densidad,d$density))+0.1), xlim=c(min(c(S[S>0])), max(c(S[S>0],x))),
#      xlab="S", main = paste("")) #main = paste("20 bars")
# lines(spline(x,densidad),col="black",lwd =3,lty=2) 
# lines(x,d12$fs,col='gray',lwd =3) #analytical distribution
# par(ps = 12,cex=1.5)
# legend(max(c(S[S>0],x[densidad>0]))-7,max(c(densidad,d$density))-0.02, c("SME","REAL"), lty=c(2,1), 
#        lwd=c(3,3),col=c("black","gray")) 
# 
# 
# ################################################################################################
# lambda=c( 36.977 -186.085  317.374  -30.322 -128.930  -77.408    6.533   72.573)
# 
# ############################################################################
# 
# 
# ########################
# #DECOMPOUNDING
# ########################
# 
# source("disentangling.R")
# source("SME_functions.R")
# 
# #####################################################################################################################
# 
# #with real distribution
# poisson_par=10
# #laplace of the gamma distribution 
# a=5
# b=15
# mu1=exp(-poisson_par*(1-(b^a)*(alpha+b)^(-a)))
# mu1
# [1] 0.022577831 0.114779965 0.220435999 0.313029725 0.389577718 0.452379160
# [7] 0.504258899 0.547589489
# 
# hu=laplace_individuals(lambda=NULL,alpha=NULL,distribution="Poisson",Lambda=poisson_par,globalZ=1e-2,mu=mu1) #0.00155
# hu
# c(0.62092132 0.78352617 0.84878521 0.88385429 0.90573081 0.92067654 0.93153345 0.93977706)
# 
# #######################################################################################################################
# #with numerical aproximation
# poisson_par=10
# #Prob of N=0
# laplace=NULL
# po=length(S[S==0])/length(S)
# for (i in 1:8)
# {
#   laplace[i]=(1/length(S))*sum(exp(-alpha[i]*S)) 
# }
# aux=(laplace-po)/(1-po)
# 
# hu=laplace_individuals(lambda=NULL,alpha=NULL,distribution="Poisson",Lambda=poisson_par,globalZ=1e-2,mu=aux) #0.00155
# hu
# hu=c(0.62363626, 0.78429340, 0.84930105, 0.88428945, 0.90611674, 0.92102470, 0.93185059, 0.94006806)
# ###########################################################################################
# #with lambda
# hu=laplace_individuals(lambda,alpha,distribution="Poisson",Lambda=poisson_par,globalZ=1e-2,mu=NULL) #0.00155
# hu
# c(0.6272 0.7875 0.8520 0.8866 0.9081 0.9228 0.9334 0.9414)
# 
# #######################################################################################################
# #u=c(0.6236, 0.7842, 0.8493, 0.8842, 0.9061, 0.9210, 0.9318, 0.9400) #with dats
# #u=c(0.6271, 0.7874, 0.8520, 0.8866, 0.9081, 0.9227, 0.9333, 0.9414) #with maxent
# #u=c(0.6209, 0.7835, 0.8487, 0.8838, 0.9057, 0.9206, 0.9315, 0.9397) # with real
# ########################################################################################################
# 
# 
# 
# #BB method...
# tol=4e-4 #9e-4() #4e-4(usando laplace muestral)    # (usando integrales)7e-4
# lambda_p=rep(5,8)
# source("opt_methods.r")
# hh=BB_modified(lambda_p,alpha,hu,M=5000,tolerance=tol,step=2,correction=0,GLL=1,
#                sigma1=0.1,sigma2=0.5,phi_max=1e20,phi_min=1e-20,N=150)
# if (h$error==000 | h$error==-999 | h$error==-666) error(h$error)
# lambda_p=hh$lambda_Opt
# 
# #lambda_p_indirect=c(48.7133399,  -1.1026692, -22.8450649, -29.3646340, -30.2645562, -29.0893322,
# #-27.1814058, -25.0788632)
# 
# #lambda_p_direct=c(37.3700503,  -6.2161901, -17.2004769, -19.9887427, -20.1397308, -19.3815749,
# #-18.3276438, -17.2144317)
# 
# #lambda_p_real=c(303.19499, -2070.21391,  2018.23140,  1567.55905,   513.00035,  -395.44164,
# #-1073.81803, -1560.76714)
# 
# #density and graphs of S and estimated
# source("SME_functions.r")
# x=seq(0,1,length=50) #x o S.
# den=SME_density(lambda_p,alpha,x,S=NULL,globalZ=0.001055) #0.001055
# Z=den$Z
# fS=den$densidadS
# plot(x,fS,type="o")
# 
# #density and graphs of S and estimated
# source("SME_functions.r")
# x=seq(0,1,length=50) #x o S.
# den=SME_density(lambda_p,alpha,x,S=NULL,globalZ=1e-1) #0.001055
# Z=den$Z
# fS_direct=den$densidadS
# #plot(x,fS,type="o")
# 
# ############################################################################################
# #density and graphs of S and estimated
# source("SME_functions.r")
# x=seq(0,1,length=50) #x o S.
# den=SME_density(lambda_p_real,alpha,x,S=NULL,globalZ=1e-1) #0.001055
# Z=den$Z
# fS_real=den$densidadS
# plot(x,fS,type="o")
# 
# 
# x11()
# #cargar data
# source("disentangling.r")
# indi=read.table(file ="individuals6.dat", header = FALSE)
# simullog=4/10*dgamma(x,shape=5,rate=15) +6/10*dgamma(x,shape=5,rate=15)
# hist(indi[,1],freq = FALSE,breaks=15,main="",ylab="density",xlab="Severity",ylim=c(0,max(simullog)))
# lines(x,simullog,ylim = c(0,max(simullog)+0.1),main="",ylab="density",xlab="Severity",col="black",type="l",lwd=3.5)
# lines(x,fS_indirect,type="l",col="gray",lty=2,lwd=3.5)
# #lines(x,fS_direct,col="gray",lty=2,lwd=3.5)
# lines(x,fS_real,type="l",col="gray",lty=1,lwd=3.5)
# legend(max(x)-.5,3, c("EquivMix", "SME", "SME (Theo. Laplace)" ), lty=c(1,3,1), lwd=c(3.5,3.5,3.5),col=c("black","gray","gray")) 
# 
# 
# source("disentangling.r")
# 
# density_simple_errors_SME(lambda_p_direct,alpha,indi,globalZ=1e-5,bins=15,ejemplo=6)
# density_simple_errors_SME2(indi,lambda_p_direct,alpha,globalZ=1e-5,ejemplo=6)
# 
# density_simple_errors_SME(lambda_p_indirect,alpha,indi,globalZ=1e-5,bins=15,ejemplo=6)
# density_simple_errors_SME2(indi,lambda_p_indirect,alpha,globalZ=1e-5,ejemplo=6)
# 
# density_simple_errors_SME(lambda_p_real,alpha,indi,globalZ=1e-3,bins=15,ejemplo=6)
# density_simple_errors_SME2(indi,lambda_p_real,alpha,globalZ=1e-5,ejemplo=6)
# 
# ########################################################################################################################################
# 
# ###################################################
# #Density reconstruction with errors in the data   #
# ###################################################
# 
# dec=1.5
# dec1=rep(dec,k)
# kk1=1:k
# alpha=dec1/kk1
# 
# 
# source("opt_methods.r")
# source("SME_functions.r")
# 
# ###########################################
# OJO: REVISAR: hay que cambiar cosas aqui
# ########
# source("disentangling.r")
# llp=boot_mu(S,k,iter=100,subsample=50,replacement=0,gamma=0.05)
# llp$mu_m
# hu_mean=llp$hu_m
# hu=llp$inter
# hu
# #########################################################################################################
# #direct
# lval       lval
# vc 0.62190247 0.62397579
# vc 0.78363061 0.78469897
# vc 0.84890713 0.84956363
# vc 0.88404679 0.88457379
# vc 0.90587818 0.90636061
# vc 0.92085922 0.92122459
# vc 0.93171709 0.93202756
# vc 0.93996779 0.94024596
# 
# hu1=c(0.62190247,0.78363061,0.84890713,0.88404679,0.90587818,0.92085922,0.93171709,0.93996779)
# hu2=c(0.61396779,0.78469897,0.84956363,0.88457379,0.90636061,0.92122459
#       ,0.93202756,0.94024596)
# hu=cbind(hu1,hu2)
# #######################################################################################################
# #indirect
# lval       lval
# vc 0.62540748 0.62719939
# vc 0.78537059 0.78748673
# vc 0.85033768 0.85206418
# vc 0.88348803 0.88479882
# vc 0.90553058 0.90633158
# vc 0.92058362 0.92108295
# vc 0.93150506 0.93187698
# vc 0.93978136 0.94004536
# 
# hu1=c(0.62540748, 0.78537059, 0.85033768, 0.88348803, 0.90553058, 0.92058362, 0.93150506, 0.93978136)
# hu2=c(0.62719939,0.78748673, 0.85206418, 0.88479882, 0.90633158, 0.92108295, 0.93187698, 0.94004536)
# hu=cbind(hu1,hu2)
# hu
# 
# 
# hu1=c(0.62710748, 0.78737059, 0.85193768, 0.88618803, 0.90793058, 0.92248362, 0.93330506, 0.94138136)
# hu2=c(0.62739939, 0.78758673, 0.85216418, 0.88669882, 0.90823158, 0.92288295, 0.93347698, 0.94154536)
# hu=cbind(hu1,hu2)
# hu
# ################################################################################
# 
# 
# #(1)
# #u=c(0.6236, 0.7842, 0.8493, 0.8842, 0.9061, 0.9210, 0.9318, 0.9400) #with dats
# #u=c(0.6271, 0.7874, 0.8520, 0.8866, 0.9081, 0.9227, 0.9333, 0.9414) #with maxent
# #u=c(0.6209, 0.7835, 0.8487, 0.8838, 0.9057, 0.9206, 0.9315, 0.9397) # with real
# #hu=cbind(u-(0.0001/2)*u,u+(0.0001/2)*u)
# #hu=cbind(u,u)
# 
# global<<-0.0001 
# lambda=rep(1,k)
# #fix always a tolerance of 1e-2 gives good results in general, in some cases
# #you can improved, but not always. Depends of the case.
# tol=4.07e-4 #1.07e-4 #8.82e-05  #8e-4
# source("opt_methods.r")
# #Initial values
# hh2=BB_modified(lambda,alpha,hu,nameFun="fun_mu_gamma.r",nameGrad="gradient_mu_gamma.r",
#                 M=10000,tolerance=tol,step=2,correction=0,GLL=1,
#                 sigma1=0.1,sigma2=0.5,phi_max=1e20,phi_min=1e-20,N=100)
# if (h$error==000 | h$error==-999 | h$error==-666) error(h$error)
# lambda_p2=hh2$lambda_Opt
# 
# #direct
# lambda_p2_direct=lambda_p2
# #lambda_p2_direct=c( 49.635702, -19.590095, -22.903152, -22.498922, -20.683808, -19.253131, -17.756677,
# #-16.454699)
# 
# #lambda_p2_indirect=c( -40.64574698,  313.62257226,  -27.33596673,   -0.67378206,   -0.98025470,
# #-0.68636640, -433.31003075, -268.30833843)
# 
# ##################
# # density
# #################
# 
# graphics.off()
# #density and graphs of S and estimated
# source("SME_functions.r")
# x=seq(0,1,length=50) #x o S.
# den=SME_density(lambda_p2_direct,alpha,x,S=NULL,globalZ=1e-7)
# Z=den$Z
# fS_direct=den$densidadS
# 
# graphics.off()
# #density and graphs of S and estimated
# source("SME_functions.r")
# x=seq(0,1,length=50) #x o S.
# den=SME_density(lambda_p2_indirect,alpha,x,S=NULL,globalZ=1e-7)
# Z=den$Z
# fS_indirect=den$densidadS
# 
# 
# 
# sum_e<-source("sum_exp.R")$value #sum of exponentials
# y<- seq(from=0,to=1,0.01)
# f_y=sum_e(lambda_p2,alpha,y)
# mus=NULL
# Z<-trapz(y,f_y)
# for (i in 1:length(lambda_p2))
# {
#   mus[i]<-trapz(y,(y^alpha[i])*(f_y/Z))
# }
# mu
# mus
# 
# x11()
# #cargar data
# source("disentangling.r")
# indi=read.table(file ="individuals6.dat", header = FALSE)
# hist(indi[,1],freq = FALSE,breaks=15,main="",ylab="density",xlab="Severity",)
# simullog=4/10*dgamma(x,shape=5,rate=15) +6/10*dgamma(x,shape=5,rate=15)
# lines(x,simullog,ylim = c(0,max(simullog)+0.1),main="",ylab="density",xlab="Severity",col="black",type="l",lwd=3.5)
# #lines(x,fS_indirect,,col="gray",lty=2,lwd=3.5)
# lines(x,fS_direct,col="gray",lty=2,lwd=3.5)
# legend(max(x)-0.2,2.7, c("EquivMix","SMEE"), lty=c(1,3), lwd=c(3.5,4),col=c("black","gray")) 
# 
# source("disentangling.r")
# density_simple_errors_SME(lambda_p2_direct,alpha,indi,globalZ=1e-5,bins=15,ejemplo=6)
# density_simple_errors_SME2(indi,lambda_p2_direct,alpha,globalZ=1e-5,ejemplo=6)
# 
# density_simple_errors_SME(lambda_p2_indirect,alpha,indi,globalZ=1e-5,bins=15,ejemplo=6)
# density_simple_errors_SME2(indi,lambda_p2_indirect,alpha,globalZ=1e-5,ejemplo=6)
# 
