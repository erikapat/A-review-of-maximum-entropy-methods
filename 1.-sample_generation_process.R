#################################################################
#DESCRIPTION:                                                   #
# The sample generarion process
# #################################################################

rm(list=ls()) #remover cualquier variable del espacio de trabajo :)
source('src/simulations.R')
require(dplyr)
require(data.table)

# S1
 simul_compound('ejemplo2_paper_2018', ele = 4, me = 6, sigma = .5, M = 2000, seed = 5)
 simul_compound('ejemplo2_paper_2018_100', ele = 4, me = 6, sigma = .5, M = 100, seed = 5)
 simul_compound('ejemplo2_paper_2018_500', ele = 4, me = 6, sigma = .5, M = 500, seed = 1)


 #S2
 simul_compound_gamma('ejemplo2_paper_2018_gamma', ele = 8, a = 350, b= 3, M = 2000, seed = 5)
 name.file <- 'ejemplo2_paper_2018_gammacompounded_dist.dat'
 S2 <- fread(paste0('data/', name.file)) %>% as.matrix() %>%  as.vector()
 hist(S2)
 
 name.file <- 'ejemplo2_paper_2018frequency.dat'
 n1 <- fread(paste0('data/', name.file)) %>% as.matrix() %>%  as.vector()
 name.file <-'ejemplo2_paper_2018_gammafrequency.dat'
 n2 <- fread(paste0('data/', name.file)) %>% as.matrix() %>%  as.vector()
 
frequency <-  c(n1, n2)
write.table(frequency, file = paste0('data/', 'frequency_two_source_risk.dat'), row.names = FALSE, col.names = FALSE)
#-------------------------------------------------------------------------------------------------
# negatives distribution
# set.seed(1)
# par(mfrow= c(1, 3))
# S = rnorm(110, -2, .2) 
# S = S*(-1)
# write.table(S, file = paste0('data/', 'negative_dist.dat'), row.names = FALSE, col.names = FALSE)
# hist(S, breaks = 10)
# min <- min(S) + .1*min(S)
# hist(S- min, breaks = 10)
# hist(S- min + min , breaks = 10)
# par(mfrow= c(1, 1))
