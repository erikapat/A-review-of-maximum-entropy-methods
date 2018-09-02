# Graphs of the papers

rm(list=ls())
source('utils/Rlibraries.R')

# 2000 points ----------------------------------------------------------------------------------------------------------------------------
b = 20 # number of bins

name.file <- 'ejemplo2_paper_2018_500compounded_dist' ###
S_new <- fread(paste0('data/', name.file, '.dat')) %>% as.matrix() %>%  as.vector()
S_new <- data.frame(S = S_new)

name.file <- 'ejemplo2_paper_2018compounded_dist' ###
#name.file <- 'ejemplo2_paper_2018_100compounded_dist' ###
S    <- fread(paste0('data/', name.file, '.dat')) %>% as.matrix() %>%  as.vector()
S <- data.frame(S = S)

output.dir <- paste0("output/", name.file, "/SME")
densi_SME  <- read.table(file =  paste0(output.dir, '/density_maxent_', nrow(S), '.txt')) %>% mutate(name = 'SME')
output.dir <- paste0("output/", name.file, "/SMEE2")
densi_SMEE  <- read.table(file =  paste0(output.dir, '/density_maxent_', nrow(S), '.txt'))  %>% mutate(name = 'SMEE')
output.dir <- paste0("output/", name.file, "/MEM_Poisson")
densi_MEM_Poisson  <- read.table(file =  paste0(output.dir, '/density_maxent_', nrow(S), '.txt'))  %>% mutate(name = 'MEM Poisson')

data = rbind(densi_SME, densi_SMEE, densi_MEM_Poisson)


num_digits <- nchar(as.character(trunc(max(S))))

if( num_digits < 4){
  data <- data %>% filter(densidad > 1e-6)
}else{
  data <- data %>% filter(densidad > 1e-8)
}

densi_gg <- ggplot(S_new) + geom_histogram(aes(x = S, y = ..density..), bins = b, color = 'lightgray', fill = 'gray40')  # 
densi_gg <- densi_gg + geom_line(data = data, aes(x = x, y = densidad, color = name, linetype = name), size = 1) 
#densi_gg <- densi_gg + geom_point(data = data, aes(x = x, y = densidad, color = name), size = .1)
#densi_gg <- densi_gg + geom_freqpoly(data = df.density, aes(x = data, y = ..density..),  bins = sep, color = 'black', linetype = 5)
densi_gg <- densi_gg +  theme_tq() 
densi_gg <- densi_gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12, face="bold"))
densi_gg <- densi_gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12, face="bold"))
densi_gg <- densi_gg + theme(legend.title = element_text(colour="blue", size=0, 
                                                         face="bold"), 
                             legend.text = element_text(colour="black", size=16)) 
densi_gg <- densi_gg + scale_color_manual(values = c( 'SME' = 'black',
                                                      'SMEE' = 'gray30',
                                                      'MEM Poisson' = 'darkgray'))
densi_gg <- densi_gg +   guides(col = guide_legend(ncol = 1)) 
densi_gg <- densi_gg + ggtitle("DENSITY DISTRIBUTION")   + xlab("S") + ylab("density")  
densi_gg <- densi_gg + theme(plot.title = element_text(size = 15, face = "bold"), legend.position=c(.8,.9))
densi_gg
dir <- paste0(getwd(), '/1.-TEXTO_PAPER/NUMERICAL_CASE_PAPER/')
if(!is.null(dir)) ggsave(paste0(dir, '/density_comparison', nrow(S), '.png'))

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

SS <- sort(S$S)
output.dir <- paste0("output/", name.file, "/SME")
F_True  <- read.table(file =  paste0(output.dir, '/F_True_', nrow(S), '.txt'))
F_SME  <- read.table(file =  paste0(output.dir, '/F_maxent_', nrow(S), '.txt')) %>% mutate(diff = F - F_True$F, name = 'SME')
output.dir <-  paste0("output/", name.file, "/SMEE2")
F_SMEE  <- read.table(file =  paste0(output.dir, '/F_maxent_', nrow(S), '.txt'))  %>% mutate(diff = F - F_True$F, name = 'SMEE')
output.dir <-  paste0("output/", name.file, "/MEM_Poisson")
F_MEM_Poisson  <- read.table(file =  paste0(output.dir, '/F_maxent_', nrow(S), '.txt'))  %>% mutate(diff = F - F_True$F, name = 'MEM Poisson')

data = rbind(F_SME, F_SMEE, F_MEM_Poisson)


gg <- ggplot(data) + geom_line(aes(x = x, y = diff, linetype = factor(name), color = factor(name)), size = 1) + ylim(c(-.11, .11)) + geom_hline(yintercept = 0)
gg <- gg + geom_hline(yintercept = -.10, color = 'red', linetype = 2) + geom_hline(yintercept = .10, color = 'red', linetype = 2)  #erase?
gg <- gg + theme_tq() 
gg <- gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12))
gg <- gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12))
gg <- gg + theme(legend.title = element_text(colour="black", size=0, face="bold"),
                 legend.text = element_text(colour="black", size=16)) 
gg <- gg + ggtitle("CALIBRATION DIAGRAM") + xlab("S") + ylab("CDF True - CDF Maxent ") + theme(plot.title = element_text(size = 15, face = "bold"))
gg <- gg + scale_color_manual(values = c( 'SME' = 'black',
                                                      'SMEE' = 'gray30',
                                                      'MEM Poisson' = 'darkgray'))
gg <- gg + theme(plot.title = element_text(size = 15, face = "bold"), legend.position=c(.8,.8))
gg

dir <- paste0(getwd(), '/1.-TEXTO_PAPER/NUMERICAL_CASE_PAPER/')
if(!is.null(dir)) ggsave(dir, '/calibration_diagram_', nrow(S), '.png'))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

source('src/ValidationFunctions.R')
# Measures
metrics(F_SME$F, F_True$F)
metrics(F_SMEE$F, F_True$F)
metrics(F_MEM_Poisson$F, F_True$F)

