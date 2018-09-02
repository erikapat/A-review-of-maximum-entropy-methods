
# FitDistributionsPrincipal

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

source('src/FitDistributionSimulations.R')



TAM = 10

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
#---------------------------------------------------------------------------
S = S1 + S2
S <- S[S > 0] 


  # simulFitBurr(S, tam = TAM, Mnumber=200, minVal=3.5, maxVal=7, interVals=30)
  # densi <- read.table(file = paste0('/home/erika/Dropbox/4ErikaIII (continuacion)/paper (5) scaling/code/FitDistributions/Burr_density_', TAM, '.dat'), header = T)
  # simGraph(densi,S,minVal=3.5,maxVal=7,porc=1.83)
  
  #densi.estimated = simulFitBurr(S, tam = 35, Mnumber=200, minVal= 0.01, maxVal= 12000, interVals= 100)
  #densi.estimated = simulFitGEV(S,tam = 11,Mnumber= 200, minVal= 0.01, maxVal= 12000, interVals= 100)

  densi.estimated = simulFitLogNormal(S,tam= 5,Mnumber=200, minVal= 0.01, maxVal= 12000, interVals=100)

  output.dir <- 'output/2nd_level/'
  densi.true <- read.table(paste0(output.dir, '/density_mat_', 2000, '.txt'))
  sampling_comparison(densi.true, densi.estimated, msn, limits = c(0, 6e-4))

