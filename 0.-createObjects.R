
rm(list=ls())
source('src/opt_methods.R')
source('utils/Rlibraries.R')


# Create objects


input            <- NULL
input$method     <- 'SMEE2' #'MEM_Poisson'  #'MEM_Exponential' #'MEM_Poisson' #'SME' 'MEM_Poisson' 'SMEE2'
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
input$alpha    <- alpha_function(input$k)
input$lambda   <-  rep(1, input$k)

saveRDS(input, file = 'input_objects/SMEE_OBJECT.RDS')

#-----------------------------------------------------------------------------------------------------------------------------------------------------
input            <- NULL
input$method     <- 'SME' #'MEM_Poisson'  #'MEM_Exponential' #'MEM_Poisson' #'SME' 'MEM_Poisson' 'SMEE2'
input$iterations <- 1000
input$tol        <- 1e-6
input$step       <- 2
input$GLL        <- 1
input$sigma[1]   <- 0.1
input$sigma[2]   <- 0.5
input$phi_max1   <- 1e20
input$phi_min1   <- 1e-20
input$seed       <- 1
input$k          <- 8 # number of momnets
input$alpha    <- alpha_function(input$k)
input$lambda   <-  rep(1, input$k)

saveRDS(input, file = 'input_objects/SME_OBJECT.RDS')

#--------------------------------------------------------------------------------------------------------------------------------------------------------

input            <- NULL
input$method     <- 'MEM_Poisson' #'MEM_Poisson'  #'MEM_Exponential' #'MEM_Poisson' #'SME' 'MEM_Poisson' 'SMEE2'
input$iterations <- 1000
input$tol        <- 1e-6
input$step       <- 2
input$GLL        <- 1
input$sigma[1]   <- 0.1
input$sigma[2]   <- 0.5
input$phi_max1   <- 1e20
input$phi_min1   <- 1e-20
input$seed       <- 1
input$k          <- 8 # number of momnets
input$alpha    <- alpha_function(input$k)
input$lambda   <-  rep(1, input$k)

saveRDS(input, file = 'input_objects/MEM_Poisson_OBJECT.RDS')