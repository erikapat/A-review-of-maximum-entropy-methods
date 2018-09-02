

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
#---------------------------------------------------------------------------
 S = S1 + S2
dt_ini <- S

#--------------------------------------------------------------------------------
#Calculation ...

trans   <- transformation(dt_ini)
product <- trans$product
dt      <- trans$dt
#-------------------------------------------------------------------------------------------
N = 1000
method = 'SME'

# output.dir <- dir_output_creation(name.file, method)
# fig_dir    <- dir_fig_creation(name.file, method)
output.dir <- 'output/2nd_level/'
source('src/visual_functions.R')
lambda           <- read.table(paste0(output.dir,  '/lambda_mat_', N, '.txt'))
densi.estimated  <- read.table(paste0(output.dir,  '/density_mat_', N, '.txt'))
#densi.true <- read.table(paste0('data/ejemplo_paper_2018true.densi.dat'))
output.dir <- 'output/2nd_level/'
densi.true <- read.table(paste0(output.dir, '/density_mat_', 2000, '.txt'))

msn = method
gg <- sampling_comparison(densi.true, densi.estimated, msn)
gg
dir <- '/home/erika/Dropbox/4ErikaIII (continuacion)/(R) CODIGO R MAXENT NEW/1.-TEXTO_PAPER/NUMERICAL_CASE_PAPER/'
if(!is.null(dir)) ggsave(paste0(dir, '/density_', method, '_', N, '.png'))

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#VaR/TVaR
gamma = .90
# VaR  <- read.table(file =  paste0(output.dir, "/VaR_SME_TAM_", N, '_gamma_', gamma,  ".txt"))
# TVaR <- read.table(file =  paste0(output.dir, "/TVaR_SME_TAM_", N, '_gamma_', gamma,  ".txt"))
source('src/visual_functions_sample_variability.R')
d_value <- draw_risk_measures_VaR_TVaR(gamma, output.dir, fig_dir)
dir <- '/home/erika/Dropbox/4ErikaIII (continuacion)/(R) CODIGO R MAXENT NEW/1.-TEXTO_PAPER/NUMERICAL_CASE_PAPER/'
d_value$d_TVaR
if(!is.null(dir)) ggsave(paste0(dir, '/box_plot_', 'TVaR', '.png'))
d_value$d_VaR
if(!is.null(dir)) ggsave(paste0(dir, '/box_plot_', 'VaR', '.png'))





