

draw_risk_measures_VaR_TVaR <- function(gamma = .95, output.dir, fig_dir){
  
  #tam   <- c(10)
  tam   <- c(10, 20, 50, 200, 500, 1000)
  #gamma <- .9
  
  VaR <- NULL
  TVaR <- NULL
  tam.x <- NULL
  for (i in 1:length(tam)){
    
    print(i)
    # load(file =  paste0(output.dir, "VaR_SME_TAM_", tam[i], '_gamma_', gamma,  ".RData"))
    # load( file = paste0(output.dir, "TVaR_SME_TAM", tam[i], '_gamma_', gamma,  ".RData"))
    m_VaR =  as.numeric(as.matrix(read.table(file =  paste0(output.dir, "/VaR_SME_TAM_", tam[i], '_gamma_', gamma,  ".txt"))))
    m_TVaR = as.numeric(as.matrix(read.table(file = paste0(output.dir, "/TVaR_SME_TAM_", tam[i], '_gamma_', gamma,  ".txt"))))
    
    VaR  <- c(VaR, m_VaR)
    TVaR <- c(TVaR, m_TVaR)
    tam.x  <- c(tam.x, rep(tam[i], length(m_VaR)))
    
    
  }
  mat_int2 <- as.data.table(cbind(tam.x, VaR, TVaR))
  setnames(mat_int2, names(mat_int2), c('tam', 'VaR', 'TVaR'))
  
  
  d <- ggplot(mat_int2, aes(x = as.factor(tam), y = VaR)) + geom_boxplot()
  d <- d + xlab('Sample size') + ylab('VaR')  + ylim(c(0, 1.1*max(mat_int2$VaR)))#+ ylim(c(min(mat_int2$g_distorted)-.01, max(mat_int2$g_distorted) + .01))
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
  d_VaR <- d
  #ggsave(d, file = paste0(fig_dir, 'VaR_','.png'), width=160, height=100, units="mm")
  
  d <- ggplot(mat_int2, aes(x = as.factor(tam), y = TVaR)) + geom_boxplot()
  d <- d + xlab('Sample size') + ylab('TVaR')  + ylim(c(0, 1.1*max(mat_int2$TVaR)))
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
  d_TVaR <- d
  
  #ggsave(d, file = paste0(fig_dir, 'TVaR_','.png'), width=160, height=100, units="mm")

  return(list(d_TVaR = d_TVaR, d_VaR = d_VaR))
}


