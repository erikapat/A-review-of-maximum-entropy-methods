

rm(list=ls())
source('utils/Rlibraries.R')
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
source("src/decompounding_functions.R")
#--------------------------------------------------------------------------------------------------------------------------------------------------
# EJEMPLO DEL PAPER ---------------------------------------------------------------------------------------------------------------------------------------
#iNDIVIDUAL LOSSES
name.comp            <- 'ejemplo2_paper_2018compounded_dist'
name.file            <- 'ejemplo2_paper_2018.individal_lognormal' ###
name.file.validation <- 'ejemplo2_paper_2018.individal_lognormal' ### 'ejemplo2_paper_2018_500.individal_lognormal' ###

#El lambda óptimo es el resultado del escalamiento, 
#por lo que hay qu aplicar la misma escala a este paso.
poisson_par <- 5

S       <- fread(paste0('data/', name.file, '.dat')) %>% as.matrix() %>%  as.vector()
dt_ini  <- S
trans   <- transformation(dt_ini)
product <- trans$product
dt      <- trans$dt

method   = 'MEM_Poisson' # SMEE2, MEM_Poisson, SME
method.i = 'MEM_Poisson' # SMEE2, MEM_Poisson, SME

output.dir  <- paste0("output/", name.comp, "/", method)
lambda      <- read.table(file =  paste0(output.dir, '/lambda_opt_', length(S), '.txt'))
output.dir  <- paste0("output/", name.comp, "/", method)
lambda      <- read.table(file =  paste0(output.dir, '/lambda_opt_', length(S), '.txt'))
output.dir  <- paste0("output/", name.comp, "/", method)
lambda      <- read.table(file =  paste0(output.dir, '/lambda_opt_', length(S), '.txt'))


# SME SMEE, MEM_PoissoN
lambda_opt = lambda %>%  as.matrix() %>% as.vector()
input <- readRDS(paste0("input_objects/", method.i, "_OBJECT.RDS")) # SMEE_OBJECT SME_OBJECT MEM_Poisson_OBJECT.RDS
input$S <- dt

#----------------------
# Decompounding
#-----------------------
hu = laplace_individuals(lambda=lambda_opt,alpha=input$alpha,distribution="Poisson", globalZ=1e-2, Lambda = poisson_par)
hu

input$mu     <- hu
lambda_opt_indi   <- maxent_types(input)

#------------------------------------------------------------------------------------------------------
# VALIDATION
# values to check...
# I have to prove the algorithm in new values

new_data <- fread(paste0('data/', name.file.validation, '.dat')) %>% as.matrix() %>%  as.vector()
new_S    <- apply_transformation_to_other_data(dt_ini, new_data)$dt
x        <- generate_x_values(new_S, num = 500)
densidad <- density.values(x, input, lambda_opt_indi,  globalZ = 1e-2)

alpha    <- input$alpha
# CDF ------------------------------------------------------------------------------------------------------------------------------
x   <- new_S
F_e <- empirical_CDF(x)$F_true #EMPIRIC DISTRIBUTION FUNCTION

if(input$method == "MEM_Poisson"){
  F_maxent <- MEMpoisson_distribution(x, lambda_opt_indi, alpha, N = 200, ita = 2, meth=4)
}else{
  F_maxent <- F_diff(x, lambda_opt_indi, alpha)$F_maxent
}
#-----------------------------------------------------------------------------------------------------------------------------------------
# Transformed data: can be useful later

S_scaled        <- dt
densidad_scaled <- densidad

#------------------------------------------------------------------------------------------
# reverse the Transformations
densidad <- inverse_transformation(dt_ini, new_S, densidad)
F_maxent = F_inverse_transformation(dt_ini, new_S, F_maxent)
F_e = F_inverse_transformation(dt_ini, new_S, F_e)


data_decompounding <- NULL
data_decompounding <- data.frame(parameter = poisson_par, RMSE = metrics(F_maxent$F, F_e$F)$RMSE)
#write.table(data_decompounding, file = paste0('output/decompounding_RMSE_', method.i, '.txt'), quote = F, append = T, col.names = T, row.names = F)
#data_decompounding <- data.frame(pos = poisson_par, RMSE = metrics(F_maxent$F, F_e$F)$RMSE)
#-vALIDATION ----------------------------------------------------------------------------------------------------------------------------------


out_dir <- paste0(name.comp, '/individuals/')
dir_output_creation(out_dir)
fig_dir <- paste0('output/', name.comp, '/individuals/', input$method, '/')
# write.table(densidad, file =  paste0(fig_dir, 'density_maxent_severity', length(dt_ini), '.txt'))
# write.table(lambda_opt, file =  paste0(fig_dir, '/lambda_opt_indi', length(dt_ini), '.txt'))
# write.table(F_maxent, file =  paste0(fig_dir, '/F_maxent_', length(dt_ini), '.txt'))
# write.table(F_e, file =  paste0(fig_dir, '/F_True_', length(dt_ini), '.txt'))
# DENSITY
g1 <- graph_histReal_curveEstimated(densidad, new_data, dir = fig_dir)
#DISTRIBUTION FUNCTION PLOT
g2 <- DistributionPlot(F_e$x, F_e$F, F_maxent$F, dir = fig_dir)
#REALIABILITY DIAGRAM
g3 <-DrawRealiabilityDiagram(F_maxent$F, F_e$F, F_e$x, dir = fig_dir)
#CALIBRATION DIAGRAM
g4 <- DrawCalibrationDiagram(F_maxent$F, F_e$F, F_e$x, dir = fig_dir)

multiplot(g1, g2, g3, g4, cols = 2)


