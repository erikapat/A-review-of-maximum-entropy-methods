#-------------------------------------------------------------------------------------------------------------

rm(list=ls()) #remover cualquier variable del espacio de trabajo :)
source('utils/Rlibraries.R')
source('utils/input_verifications.R')
source('src/actuarial_functions.R')
source('src/statistical_functions.R')
source('src/maxent_algorithm.R')
source('src/SME_functions.R')


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

method = 'SME' 
input <- readRDS(paste0('input_objects/', method, '_OBJECT.RDS')) # SMEE_OBJECT SME_OBJECT MEM_Poisson_OBJECT.RDS

alpha <- c(1.5000000, 0.7500000, 0.5000000, 0.3750000, 0.3000000,
           0.2500000, 0.2142857, 0.1875000)

tam   <- c(10, 20, 50, 200, 500, 1000)
gamma   <- .9
a       <- 1
factor = 1000
len = 200

mat_int <- NULL
output.dir <- 'output/2nd_level/'

#extreme limits
min_limit_int <- 0
max_limit_int <- (1.8*max(S))/factor #to escale the maximum value possible (to calculate densities with lambda)

i <- 1
for (i in 1:length(tam)){
  
  print(tam[i])
  N     <- tam[i]
  mat   <- read.table(paste0(output.dir,  '/lambda_mat_', N, '.txt'))
  densi <- read.table(paste0(output.dir, '/density_mat_', N, '.txt'))
  mat   <- as.matrix(mat)
  
  sdx  <- standard_dev(mat, densi)
  mea  <- expectated_value(mat, densi)
  
  #stop-loss with a cap----------------------------------------------------------------------------------------------------------------------------
  #Formula (3) en el paper: Sample dependence of risk premia
  
  #stop-loss with a cap
  # load(paste0("input/VaR_SME_TAM_", tam[i], "_gamma_", gamma, ".RData"))
  # load(paste0("input/TVaR_SME_TAM", tam[i], "_gamma_", gamma, ".RData"))
  m_VaR <- read.table(file =  paste0(output.dir, "/VaR_SME_TAM_", N, '_gamma_', gamma,  ".txt"))
  m_TVaR <- read.table(file =  paste0(output.dir, "/TVaR_SME_TAM_", N, '_gamma_', gamma,  ".txt"))
  
  d    <- as.vector(as.matrix(m_VaR))
  K    <- as.vector(as.matrix(m_TVaR))
  mat_t <- t(mat)
  integral_cap <- cap_function(mat_t, input, factor, d, K, len, max_limit_int)
  
  
  # risk premia (1): root ------------------------------------------------------------------------------------------------------------------------------------
  integral_premia <- risk_premium(densi, a, type = 'root', factor = factor)
  #boxplot(integral)
  
  # risk premia (2): exp ------------------------------------------------------------------------------------------------------------------------------------
  
  integral_premia2 <- risk_premium(densi, type = 'exp', factor = factor)
  
  #g-distorted risk price
  F_ <- read.table(file =  paste0(output.dir, "/CDF_", N, ".txt"))
  val_g            <- g_distorted_risk_price(densi, F_, tam[i], a)
  g_distorted      <- as.vector(val_g$g_distorted)
  #g_distorted_sqrt <- as.vector(val_g$g_distorted_sqrt)
  
  mat <- cbind(rep(tam[i], length(integral_cap)),  integral_cap, integral_premia, integral_premia2, g_distorted) #, g_distorted_sqrt
  
  mat_int <- rbind(mat_int, mat)
  
}  
#colMeans(mat_int[M== 1000,], na.rm = T)
#-----------------------------------------------------------------------------------------------------------------------------------------------------
mat_int        <- as.data.table(mat_int)
names(mat_int) <- c('M', 'integral_cap', 'intregral_premia', 'integral_premia2',  'g_distorted') #,  'g_distorted_sqrt')

#cap
mat_int2 <- mat_int #[integral_cap < .5, ]
d <- ggplot(mat_int2, aes(x = as.factor(M), y = integral_cap)) + geom_boxplot()
d <- d + xlab('Sample size') + ylab('stop-loss with a cap') 
d <- d + theme(axis.title.y = element_text(size = rel(1.5)))
d <- d + theme(axis.text.y = element_text(size = rel(1.8)))
d <- d + theme(axis.text.x = element_text(size = rel(1.8)))
d <- d + theme(axis.title.x = element_text(size = rel(1.5)))
d <- d + theme(legend.text = element_text(size = 15))
d <- d + theme(legend.title = element_text(size = 15))
d <- d +  theme(
  panel.background = element_rect(fill="white") ,
  panel.grid.minor.y = element_line(size=3),
  panel.grid.major = element_line(colour = "lightgray"),
  plot.background = element_rect(fill="white")
)
d
dir <- '/home/erika/Dropbox/4ErikaIII (continuacion)/(R) CODIGO R MAXENT NEW/1.-TEXTO_PAPER/NUMERICAL_CASE_PAPER/'
if(!is.null(dir)) ggsave(paste0(dir, '/stop_loss_with_a_cap', '.png'))
#premium 1
#mat_int2 <- mat_int[intregral_premia < .5, ]
mat_int2 <- mat_int
d <- ggplot(mat_int2, aes(x = as.factor(M), y = intregral_premia)) + geom_boxplot()
d <- d + xlab('Sample size') + ylab('risk premia')  #+ ylim(c(3.5, 5.5))
#d <- d + scale_y_continuous(limits = c(3.5, 5))
d <- d + theme(axis.title.y = element_text(size = rel(1.5)))
d <- d + theme(axis.text.y = element_text(size = rel(1.8)))
d <- d + theme(axis.text.x = element_text(size = rel(1.8)))
d <- d + theme(axis.title.x = element_text(size = rel(1.5)))
d <- d + theme(legend.text = element_text(size = 15))
d <- d + theme(legend.title = element_text(size = 15))
d <- d +  theme(
  panel.background = element_rect(fill="white") ,
  panel.grid.minor.y = element_line(size=3),
  panel.grid.major = element_line(colour = "lightgray"),
  plot.background = element_rect(fill="white")
)
d
if(!is.null(dir)) ggsave(paste0(dir, '/risk_premia_Integral', '.png'))
#premium 2
#mat_int2 <- mat_int[intregral_premia < .5, ]
mat_int2 <- mat_int
d <- ggplot(mat_int2, aes(x = as.factor(M), y = integral_premia2)) + geom_boxplot()
d <- d + xlab('Sample size') + ylab('risk premia')  #+ ylim(c(3.5, 5.5))
d <- d + theme(axis.title.y = element_text(size = rel(1.5)))
d <- d + theme(axis.text.y = element_text(size = rel(1.8)))
d <- d + theme(axis.text.x = element_text(size = rel(1.8)))
d <- d + theme(axis.title.x = element_text(size = rel(1.5)))
d <- d + theme(legend.text = element_text(size = 15))
d <- d + theme(legend.title = element_text(size = 15))
d <- d +  theme(
  panel.background = element_rect(fill="white") ,
  panel.grid.minor.y = element_line(size=3),
  panel.grid.major = element_line(colour = "lightgray"),
  plot.background = element_rect(fill="white")
)
d
if(!is.null(dir)) ggsave(paste0(dir, '/risk_premia_2', '.png'))

mat_int2 <- mat_int
d <- ggplot(mat_int2, aes(x = as.factor(M), y = g_distorted)) + geom_boxplot()
d <- d + xlab('Sample size') + ylab('g- distorted risk price')  #+ ylim(c(3.5, 5))
d <- d + theme(axis.title.y = element_text(size = rel(1.5)))
d <- d + theme(axis.text.y = element_text(size = rel(1.8)))
d <- d + theme(axis.text.x = element_text(size = rel(1.8)))
d <- d + theme(axis.title.x = element_text(size = rel(1.5)))
d <- d + theme(legend.text = element_text(size = 15))
d <- d + theme(legend.title = element_text(size = 15))
d <- d +  theme(
  panel.background = element_rect(fill="white") ,
  panel.grid.minor.y = element_line(size=3),
  panel.grid.major = element_line(colour = "lightgray"),
  plot.background = element_rect(fill="white")
)
d
if(!is.null(dir)) ggsave(paste0(dir, '/g_distorted', '.png'))


# mat_int2 <- mat_int
# d <- ggplot(mat_int2, aes(x = as.factor(M), y = g_distorted_sqrt)) + geom_boxplot()
# d <- d + xlab('Sample size') + ylab('g- distorted risk price')  + ylim(c(min(mat_int2$g_distorted_sqrt)-.01, max(mat_int2$g_distorted_sqrt) + .01))
# d <- d + theme(axis.title.y = element_text(size = rel(1.5)))
# d <- d + theme(axis.text.y = element_text(size = rel(1.8)))
# d <- d + theme(axis.text.x = element_text(size = rel(1.8)))
# d <- d + theme(axis.title.x = element_text(size = rel(1.5)))
# d <- d + theme(legend.text = element_text(size = 15))
# d <- d + theme(legend.title = element_text(size = 15))
# d <- d +  theme(
#   panel.background = element_rect(fill="white") ,
#   panel.grid.minor.y = element_line(size=3),
#   panel.grid.major = element_line(colour = "lightgray"),
#   plot.background = element_rect(fill="white")
# )
# d
# 
# ggsave(d, file = paste0('figures_2017/g_distorted_sqrt', max_limit_int,'.png'), width=160, height=100, units="mm")


