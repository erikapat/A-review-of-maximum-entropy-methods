

rm(list=ls())

method.i = 'SMEE'
dt.SMEE <- read.table(file = paste0('output/decompounding_RMSE_', method.i, '.txt'), header = T) %>%  mutate(method = 'SMEE', parameter = as.numeric(parameter), 
                                                                                                             RMSE = as.numeric(RMSE)) 

method.i = 'SME'
dt.SME <- read.table(file = paste0('output/decompounding_RMSE_', method.i, '.txt'), header = T) %>%  mutate(method = 'SME', parameter = as.numeric(parameter), 
                                                                                                            RMSE = as.numeric(RMSE))

method.i = 'MEM_Poisson'
dt.MEM <- read.table(file = paste0('output/decompounding_RMSE_', method.i, '.txt'), header = T) %>%  mutate(method = 'MEM_Poisson', parameter = as.numeric(parameter), 
                                                                                                            RMSE = as.numeric(RMSE))

data <- rbind(dt.SMEE, dt.SME, dt.MEM)
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
gg <- ggplot(data) + geom_line(aes(x = parameter, y = RMSE, linetype = factor(method), color = factor(method)), size = 1) 
gg <- gg + theme_tq() 
gg <- gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12))
gg <- gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12))
gg <- gg + theme(legend.title = element_text(colour="black", size=0, face="bold"),
                 legend.text = element_text(colour="black", size=16)) 
gg <- gg + ggtitle("") + xlab("parameter") + ylab("RMSE") + theme(plot.title = element_text(size = 15, face = "bold"))
gg <- gg + scale_color_manual(values = c( 'SME' = 'black',
                                          'SMEE' = 'gray30',
                                          'MEM_Poisson' = 'darkgray'))
gg <- gg + theme(plot.title = element_text(size = 15, face = "bold"), legend.position=c(.8,.3))
gg

dir <- paste0(getwd(), '1.-TEXTO_PAPER/NUMERICAL_CASE_PAPER/')
if(!is.null(dir)) ggsave(paste0(dir, '/parameter_poisson', '.png'))
