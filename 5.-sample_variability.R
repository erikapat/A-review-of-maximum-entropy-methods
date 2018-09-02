

rm(list=ls())
source('utils/Rlibraries.R')
source('utils/Multiple_ggplot_function.R')
source('utils/dir_functions.R')
source('utils/input_verifications.R')
#--
source('src/maxent_algorithm.R')
source('src/opt_methods.R')
source('src/SME_functions.R')
source('src/MEM_functions.R')
source('src/SMEE_functions.R')
source('src/ValidationFunctions.R')
source('src/visual_functions.R')

#ejemplo_paper_2018.individal_lognormal.dat
#ejemplo_paper_2018compounded_dist.dat
#name.file <- 'ejemplo_paper_2018compounded_dist.dat'
name.file <- 'ejemplo2_paper_2018compounded_dist.dat'
#----------------------------------------------------------------------
#name.file <- 'negative_dist.dat'
S1 <- fread(paste0('data/', name.file)) %>% as.matrix() %>%  as.vector()
S1 <- input_basic_verifications(S1)
#---------------------------------------------------------------------------
S1 <- S1[S1 > 0] # in our context zero is not important...
name.file <- 'ejemplo2_paper_2018_gammacompounded_dist.dat'
S2 <- fread(paste0('data/', name.file)) %>% as.matrix() %>%  as.vector()
S2 <- input_basic_verifications(S2)
S = S1 + S2
dt_ini  <- S[S> 0]
#---------------------
#-----------------------------------------------------------------------------
trans <- transformation(dt_ini)
product <- trans$product
dt <- trans$dt

dt %>% summary()
hist(dt, freq = F)

save.figures     <- T
input            <- NULL
input$method     <- 'SME' #'MEM_Exponential' #'MEM_Poisson' #'SME'
input$iterations <- 1000
input$tol        <- 1e-12
input$step       <- 2
input$GLL        <- 1
input$sigma[1]   <- 0.1
input$sigma[2]   <- 0.5
input$phi_max1   <- 1e20
input$phi_min1   <- 1e-20
input$seed       <- 1
input$k          <- 8 # number of momnets
input$alpha  <- alpha_function(input$k)
input$S <- dt

output.dir <- 'output/2nd_level/'
#first level should be created
if (!dir.exists(output.dir)){
  cat('[Info]: Creating directory ', output.dir, '\n')
  dir.create(output.dir)
}

#-------------------------------------------------------------------------------------------------------------------------------------------
M = 1 #NUMBER OF SAMPLES
N = 2000  # SIZE OF EACH SAMPLE (CAMBIAR) #200/10
# matrix with lambda values

registerDoParallel(cores=2)
lambda    <- foreach(i = 1:M, .combine = cbind) %dopar%  {
  set.seed(i)
  ind     <- sample(1:length(dt), N)
  input$S <- dt[ind]
  mp      <- maxent_parameters(input$S, input$k)   #aqui se deben incluir los ceros (OJO)
  #input$alpha <- mp$alpha
  input$mu     <- mp$mu
  input$lambda <- rep(1, input$k)
  lamb_opt     <- maxent_types(input)
  return(lamb_opt)
}
#stopCluster()
if (M == 1){
  lambda <- lambda %>%  as.data.frame()
}
write.table(lambda, file =  paste0(output.dir, '/lambda_mat_', N, '.txt'))

#lambda <- read.table(file =  paste0(output.dir, '/lambda_mat_', N, '.txt'))
input$S <- dt
# matrix with lambda values
alpha  <- alpha_function(input$k)
x <- generate_x_values(input$S, factor = 1.5)
mat <- foreach(i = 1:M, .combine = cbind) %do%  {
  set.seed(i)
  densidad = density.values(x, input, lambda[,i])
  if (i != 1){
    output  <- inverse_transformation(dt_ini, dt, densidad)[, 2]
  }else{
    output   = inverse_transformation(dt_ini, dt, densidad)
  }
  return(output)
}

write.table(mat, file =  paste0(output.dir, '/density_mat_', N, '.txt'))
#mat <- read.table(file =  paste0(output.dir, '/density_mat_', N, '.txt'))

#-----------------------------------------------------------------------------
# N = 1000
# M = 200
# lambda <- read.table(file =  paste0(output.dir, '/lambda_mat_', N, '.txt'))
# library(doParallel)
registerDoParallel(cores=2)
F_maxent <- foreach(i = 1:M, .combine = cbind) %dopar%  {
  if(input$method == "MEM_Poisson"){
    F_maxent    <- MEMpoisson_distribution(x, lambda[, i], input$alpha, N = 200, ita = 2, meth=4)
    F_maxent <- F_maxent[, 2]
  }else{
    F_maxent <- F_diff(x, lambda[, i], input$alpha)$F_maxent
    F_maxent <- F_maxent[, 2]
 }
  return(F_maxent)
}
F_maxent <- data.frame(x = x*product, F_maxent )
write.table(F_maxent, file =  paste0(output.dir, '/CDF_', N, '.txt'))

N = 1000
M = 200
lambda <- read.table(file =  paste0(output.dir, '/lambda_mat_', N, '.txt'))
#Calculate VaR and TVaR
alpha  <- alpha_function(input$k)
x <- generate_x_values(input$S)
source('src/VaR_TVaR_functions.R')
gamma = .90
library(doParallel)
registerDoParallel(cores=2)
VaR_TVaR_mat <- foreach(i = 1:M, .combine = cbind) %dopar%  {
  set.seed(i)
  RISK <- VaR_function_for_maxent(x, gamma, lambda[, i], alpha, globalZ=0.001)
  return(RISK)
}
VaR <- as.numeric(as.matrix(as.data.frame(VaR_TVaR_mat)[1,]))*product
TVaR <- as.numeric(as.matrix(as.data.frame(VaR_TVaR_mat)[2,]))*product
#stopCluster()

write.table(VaR, file =  paste0(output.dir, "/VaR_SME_TAM_", N, '_gamma_', gamma,  ".txt"), row.names = F, quote = F, col.names = F)
write.table(TVaR, file =  paste0(output.dir, "/TVaR_SME_TAM_", N, '_gamma_', gamma,  ".txt"), row.names = F, quote = F, col.names  = F)

source('src/visual_functions_sample_variability.R')
gg <- draw_risk_measures_VaR_TVaR(gamma = .9, output.dir, fig_dir)
gg
