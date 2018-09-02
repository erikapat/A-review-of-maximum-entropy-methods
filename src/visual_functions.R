


# VALIDATION GRAPHS --------------------------------------------------------------------------------------------------------------------------------------------


# Estimated density distribution: real (histogram) vs. maxent
# dataset: it is a data frame that contains all the information that we need
graph_histReal_curveEstimated <- function(densidad, S, dir = NULL, b = 20){
  
  S <- data.frame(S = S)
  df.density <- data.frame(x = densidad$x, densidad = densidad$densidad) 
  
  num_digits <- nchar(as.character(trunc(max(S))))
  
  if( num_digits < 4){
    df.density <- df.density %>% dplyr::filter(densidad > 1e-6)
  }else{
    df.density <- df.density %>% dplyr::filter(densidad > 1e-16)
  }
  densi_gg <- ggplot(S) + geom_histogram(aes(x = S, y = ..density..), bins = b, color = 'gray', fill = 'darkgray')  # 
  densi_gg <- densi_gg + geom_line(data = df.density, aes(x = x, y = densidad), size = .5)
  #densi_gg <- densi_gg + geom_freqpoly(data = df.density, aes(x = data, y = ..density..),  bins = sep, color = 'black', linetype = 5)
  densi_gg <- densi_gg +  theme_tq() 
  densi_gg <- densi_gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12, face="bold"))
  densi_gg <- densi_gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12, face="bold"))
  densi_gg <- densi_gg + theme(legend.title = element_text(colour="blue", size=10, 
                                                           face="bold")) 
  densi_gg <- densi_gg +   guides(col = guide_legend(ncol = 1)) 
  densi_gg <- densi_gg + ggtitle("DENSITY DISTRIBUTION")   + xlab("S") + ylab("density")  + theme(plot.title = element_text(size = 15, face = "bold"))
  
  if(!is.null(dir)) ggsave(paste0(dir, '/density.png'))

  return(densi_gg)
  
}


# Cumualtive distribution functions: TRUE vs. maxent
DistributionPlot <- function(S, F_e, F_maxent, dir = NULL){
  
  SS <- sort(S)
  df <- data.frame(S = SS, CDF.true = F_e, CDF.maxent = F_maxent)
  
  df_melt <- melt(df, id.vars = 'S')
  gg <- ggplot(df_melt, aes( x = S, y = value, color = variable)) + geom_line(aes(linetype = variable), size = 1.5)
  gg <- gg + theme_tq() 
  gg <- gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12, face="bold"))
  gg <- gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12, face="bold"))
  gg <- gg + theme(legend.title = element_text(colour="black", size=20, face="bold")) 
  gg <- gg + scale_color_manual(name="", values = c("CDF.true" = "black", "CDF.maxent" = "red"))
  #gg <- gg + ggtitle("Cumulative Distribution Function (CDF)") + xlab("S") + ylab("F(S)") 
  gg <- gg + scale_linetype_manual(name="", values = c("CDF.true" = 1, "CDF.maxent" = 2))
  gg <- gg + theme(legend.position = "top", legend.direction = "horizontal", legend.text=element_text(size= 12)) +  xlab("S") + ylab("CDF")
  
  if(!is.null(dir)){
    ggsave(paste0(dir, '/cdf.png'))
  }

  
  # VAR - TVAR - OPCIONES 
  #cat(rep("\n", 50)) #clean screen
  #diffF <-  abs(F_e - F) 
  #cat("Max. difference between distributions:", max(diffF), "\n")
  #abline(v = SS[which(diffF ==max(diffF))], col = "green")
  #abline(h = F_e[which(diffF ==max(diffF))], col = "green")
  #abline(h = F[which(diffF ==max(diffF))], col = "green")
  #text(max(SS[which(diffF ==max(diffF))]), y = mean(max(F_e[which(diffF == max(diffF))]) + max( F[which(diffF == max(diffF))])),"max")
  #cat("Max. difference between distributions:", SS[which(diffF ==max(diffF))], "\n")
  
  return(gg)
}

#REALIABILITY DIAGRAM
DrawRealiabilityDiagram <- function(F_maxent, F_e, S, dir = NULL){
  
  SS <- sort(S)
  #here we  consider only quantiles 
  qq <- round(seq(1, length(SS), length = (length(SS)/2)))
  
  if (length(qq) < 100){
    lol <- length(qq)
  }else{
    lol <- 100
  }
  
  quantiles <- sample(qq, lol, replace = FALSE, prob = NULL)
  sampF     <- F_maxent[quantiles]
  sampF_e   <- F_e[quantiles]
  df <- data.frame(quantile = quantiles, sampF_e = sampF_e, sampF = sampF)
  gg <- ggplot(df) + geom_abline(intercept = 0, color = 'gray', size = 1, linetype = 1) + geom_point(aes(x = sampF_e, y = sampF), size = 2, alpha = .6) 
  gg <- gg + theme_tq() 
  gg <- gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12, face="bold"))
  gg <- gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12, face="bold"))
  gg <- gg + theme(legend.title = element_text(colour="black", size=20, face="bold")) 
  gg <- gg + ggtitle("REALIABILITY DIAGRAM") + xlab("True Probability") + ylab("Maxent Probability") + theme(plot.title = element_text(size = 15, face = "bold"))
  gg
  
  if(!is.null(dir)){
    ggsave(paste0(dir, '/DrawRealiabilityDiagram.png'))
  }
  
  return(gg)
}


#CALIBRATION DIAGRAM
DrawCalibrationDiagram <- function(F_maxent, F_e, S, dir = NULL){
  
  SS <- sort(S)
  df <- data.frame(S = SS, dif = F_e - F_maxent)
  gg <- ggplot(df) + geom_line(aes(x = S, y = dif), size = 1, alpha = .6) + ylim(c(-.11, .11)) + geom_hline(yintercept = 0)
  gg <- gg + geom_hline(yintercept = -.10, color = 'red', linetype = 2) + geom_hline(yintercept = .10, color = 'red', linetype = 2)  #erase?
  gg <- gg + theme_tq() 
  gg <- gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12))
  gg <- gg + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12))
  gg <- gg + theme(legend.title = element_text(colour="black", size=20, face="bold")) 
  gg <- gg + ggtitle("CALIBRATION DIAGRAM") + xlab("S") + ylab("CDF True - CDF Maxent ") + theme(plot.title = element_text(size = 15, face = "bold"))
  
  cat(rep("\n", 1)) #clean screen
  cat("Max. difference between distributions:",max(abs(df$dif)),"\n")
  
  if(!is.null(dir)){
    ggsave(paste0(dir, '/CalibrationDiagram.png'))
  }
  
  return(gg)
}




#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# SAMPING DEPENDENCY GRAPHS


###############################################################################


#
# This script contains the functions used to graph multiple densities generated for methods like convolutions, maxent and others, 
#using sampling.

require(data.table)
library(plotly)
require(ggplot2)


# densi.true = densi.true
# densi.estimated = t(densi)
# msn = ''
# GRAPH
sampling_comparison <- function(densi.true, densi.estimated, msn, limits = c(0, 6e-4)){
  
  X_estimated     <- densi.estimated[, 1]
  X_TRUE          <- densi.true[, 1]
  dt.densi.true   <- as.data.table(densi.true)
  densi.estimated <- densi.estimated[, -1]
  
  #discard those row that are all zero

  #densi.estimated <- densi.estimated[which(apply(densi.estimated, 1, sum) != 0), ] #row
  #densi.estimated <- densi.estimated[which(apply(densi.estimated, 2, sum) != 0), ] #column
  #average of the reconstructions
  average  <- apply(densi.estimated, 1, mean) 
  # densities that serves like frontiers
  lin_inf <- apply(densi.estimated, 1, min) 
  lin_sup <- apply(densi.estimated, 1, max) 
  
  dt.cross <- as.data.table(cbind(as.vector(as.numeric(X_estimated)), lin_inf, lin_sup, average))
  setnames(dt.cross, names(dt.cross), c('x_est', 'lin_inf', 'lin_sup', 'average'))
  dt.cross <- interpolate_matrix(dt.cross, n = 100)
  dt.cross <- dt.cross[, col := 'True Density'][, col2 := 'Reconstructions']
  
  names(dt.densi.true) <- c('x', 'true_density')
  if(nrow(dt.densi.true) < nrow(dt.cross)){
    dt.densi.true <- interpolate_values(dt.densi.true, nrow(dt.cross))
  }
  dt.densi.true        <- dt.densi.true[, col3 := 'Average']
  
  gg <- gg_graph(dt.cross, dt.densi.true, msn, limits)
  
  return(gg)
}

#---------------------------------------------------------------------------------------------------------------

gg_graph <- function(dt.cross, dt.densi.true, msn, limits = c(0, 6e-4)){
  
  dt <- cbind(dt.cross, dt.densi.true)
  gg <- ggplot(dt, aes(x = x_est, ymin = lin_inf, ymax = lin_sup, col = col2)) +
    geom_ribbon(fill = 'gray', alpha= 0.8) 
  gg <- gg + geom_line(aes(x = x_est, y = average, col = col3), size = 1)
  gg <- gg +  geom_line(aes(x = x, y = true_density, color = col), size = 1)
  
  #gg <- gg + scale_color_manual("Click on the colors",values=c("red", "gray",  "black"))
  gg <- gg + scale_color_manual("",values=c("gray40", "gray",  "black"))
  gg <- gg + fte_theme() + theme(legend.position=  c(.9, .9))
  gg <- gg + xlab('x') + ylab('density') #+ ylim(limits)
  gg <- gg + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.title.x = element_text(color="black", size= 12, face="bold"))
  gg <- gg  + theme(axis.text.y = element_text(hjust = 1, size = 12), axis.title.y = element_text(color="black", size= 12, face="bold"))
  #gg <- gg + scale_x_continuous(labels = scales::comma)
  # gg <- ggplotly(gg, tooltip = c("y"), dynamicTicks = T)
  # gg <- ggplotly(gg, session="knitr") 
  gg
  # chart_link = api_create(gg, filename="ggplot-user-guide/2")
  # chart_link
  return(gg)
}

interpolate_matrix <- function(dt, n){
  s.n <- names(dt)
  dt.aux <- copy(dt)
  dt.n <- NULL
  if (ncol(dt) > 2){
    vect <- c('x', paste0('y', 2:ncol(dt)))
    for (i in 2:ncol(dt)){
      dt.aux <- dt[, c(1, i), with = F]
      dt.aux <- interpolate_values(dt.aux, n)
      if (i ==2){
        dt.n <- dt.aux[, c(1,2), with = F]
      }else{
        dt.n <- cbind(dt.n, dt.aux[, c(2), with = F])
      }
    }
  }
  return(dt.n)
}

interpolate_values <- function(dt, n){
  s.n       <- names(dt)
  names(dt) <- c('x', 'y')
  dt        <- spline(dt$x, dt$y, n)
  dt        <- as.data.table(cbind(dt$x, dt$y))
  names(dt) <- s.n
  
  return(dt)
}

#--------------------------------------------------------------------------------------------------------------------------------------
# alpha       <- c(1.5000000, 0.7500000, 0.5000000, 0.3750000, 0.3000000,
#                  0.2500000, 0.2142857, 0.1875000)
# lambda_5000 <- c(16240.12132,  -748.02989, -5997.50028,  1105.38763,  
#                  4745.65939,  3727.93062, -149.98752, -5356.22172)
# 
# 
# maxent_density <- function(maxent_results, ){
# 
#   x            <- sort(S[S>0])
#   F_SME        <- F_diff_new(x, lambda_5000, alpha, globalZ = 1e-5)
#   den_5000     <- F_SME$F_emp
# }

#----------------------------------------------------------------------------------------------------------------------------------------

fte_theme <- function() {
  
  require('RColorBrewer')
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = palette[9]
  color.grid.major = palette[3]
  color.axis.text = palette[6]
  color.axis.title = palette[7]
  color.title = palette[9]
  
  # Begin construction of chart
  theme_bw(base_size=9) +
    
    # Set the entire chart region to a light gray color
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    
    # Format the grid
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # Format the legend, but hide by default
    #theme(legend.position=c(0.5, 0.5)) +
    theme(legend.position= 'none')  + 
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=18,color=color.axis.title)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=20, vjust=1.25)) +
    theme(axis.text.x=element_text(size= 16,color=color.axis.text)) +
    theme(axis.text.y=element_text(size= 14,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=20,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=16,color=color.axis.title, vjust=0)) # +
  
  # Plot margins
  theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}

#-----------------------------------------------------------------------------------------------------------------------------------


