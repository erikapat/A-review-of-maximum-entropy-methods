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


#--------------------------------------------------------------
name.file <- 'ejemplo2_paper_2018compounded_dist.dat' ###
#name.file <- 'ejemplo2_paper_2018_100compounded_dist.dat' ###
dt_ini <- fread(paste0('data/', name.file)) %>% as.matrix() %>%  as.vector()
dt_ini <- input_basic_verifications(dt_ini)
#---------------------------------------------------------------------------
dt_ini <- dt_ini[dt_ini > 0] # in our context zero is not important...
#-----------------------------------------------------------------------------
trans <- transformation(dt_ini)
product <- trans$product
dt <- trans$dt

save.figures     <- T

method = 'SMEE' 
input <- readRDS(paste0('input_objects/', method, '_OBJECT.RDS')) # SMEE_OBJECT SME_OBJECT MEM_Poisson_OBJECT.RDS
input$S = dt

method = 'SMEE2' #'MEM_Poisson' #'SME' 
set.seed(input$seed)
# create folders -------------------------------------------------------------------

fig_dir     <- dir_fig_creation(name.file, method)
output.dir  <- dir_output_creation(name.file, method)
#------------------------------------------------------------------------------------------------------

# Input treatment
set.seed(input$seed)
mp    <- maxent_parameters(dt, input$k)   #aqui se deben incluir los ceros (OJO)
#input$alpha <- mp$alpha
input$mu    <- mp$mu
input$lambda <- rep(1, input$k)

#------------------------------------------------------------------------------------------------------
# Maxent algorithm: the lambdas that are necessary to regenarate the density... 
lambda_opt   <- maxent_types(input)
write.table(lambda_opt, file =  paste0(output.dir, '/lambda_opt_', length(dt_ini), '.txt'))

#------------------------------------------------------------------------------------------------------
# VALIDATION
# values to check...
# I have to prove the algorithm in new values

name.file <- 'ejemplo2_paper_2018_500compounded_dist.dat' ###
#name.file <- 'ejemplo2_paper_2018_50compounded_dist.dat' ###
new_data <- fread(paste0('data/', name.file)) %>% as.matrix() %>%  as.vector()
new_S    <- apply_transformation_to_other_data(dt_ini, new_data)$dt
x        <- generate_x_values(new_S, num = 500)
densidad <- density.values(x, input, lambda_opt)

alpha    <- input$alpha
# CDF ------------------------------------------------------------------------------------------------------------------------------
x = new_S
F_e       <- empirical_CDF(x)$F_true #EMPIRIC DISTRIBUTION FUNCTION

if(input$method == "MEM_Poisson"){
  F_maxent    <- MEMpoisson_distribution(x, lambda_opt, alpha, N = 200, ita = 2, meth=4)
}else{
  F_maxent    <- F_diff(x, lambda_opt, alpha)$F_maxent
}
#-----------------------------------------------------------------------------------------------------------------------------------------
# Transformed data: can be useful later

# S_scaled        <- dt
# densidad_scaled <- densidad

#------------------------------------------------------------------------------------------
# reverse the Transformations
densidad <- inverse_transformation(dt_ini, new_S, densidad)
F_maxent = F_inverse_transformation(dt_ini, new_S, F_maxent)
F_e = F_inverse_transformation(dt_ini, new_S, F_e)
#-vALIDATION ----------------------------------------------------------------------------------------------------------------------------------

# DENSITY
g1 <- graph_histReal_curveEstimated(densidad, new_data, dir = fig_dir)
#DISTRIBUTION FUNCTION PLOT
g2 <- DistributionPlot(F_e$x, F_e$F, F_maxent$F, dir = fig_dir)
#REALIABILITY DIAGRAM
g3 <-DrawRealiabilityDiagram(F_maxent$F, F_e$F, F_e$x, dir = fig_dir)
#CALIBRATION DIAGRAM
g4 <- DrawCalibrationDiagram(F_maxent$F, F_e$F, F_e$x, dir = fig_dir)

multiplot(g1, g2, g3, g4, cols = 2)

#-------------------------------------------------------------------------------------------------------------------------------------------

# validation

metrics(F_maxent$F, F_e$F)

#-------------------------------------------------------------------------------------------------------------------------------------------

# SAVE RESULTS
write.table(densidad, file =  paste0(output.dir, '/density_maxent_', length(dt_ini), '.txt'))

write.table(F_maxent, file =  paste0(output.dir, '/F_maxent_', length(dt_ini), '.txt'))
write.table(F_e, file =  paste0(output.dir, '/F_True_', length(dt_ini), '.txt'))
#----------------------------------------------------------------------------------------------------------------------------------------------
# Var & Tvar values
source('src/VaR_TVaR_functions.R')
#empirical values
gamma = .95
VaR_function(dt_ini, gamma, num_rep = 100, sup1=.80, inf1=0.50, seed = 10)
VaR_con_function_dist(dt_ini,gamma)
#VaR_function_table(new_data,gamma,sup1=1,inf1=.90) eliminate, this is not working
#VaR_con_function(new_data,gamma)
TVaR_function(dt_ini,gamma,sup1=1,inf1=0.80)

dt.aux <- c(0, 18)
x <- generate_x_values(dt, num = 3000)
RISK <- VaR_function_for_maxent(x,gamma,lambda_opt,alpha,globalZ=0.001)
product
RISK$VaR*product
RISK$TVaR*product
#-------------------------------------------------------------------------------------------------------------------------------------------



