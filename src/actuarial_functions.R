
#input
#d & K constants defined in the paper
# len: number of values of f(x), i.e x
# max_limit_int value that is greater enough to be infinite (superior limit in the integration)
cap_function <- function(mat,  input, factor, d, K, len, max_limit_int){
  
  alpha <- input$alpha
  densi <- NULL

  densi1 <- t(sapply(1:nrow(mat), function(i){ Maxent_Density(mat[i,], input, factor, minVal= d[i], maxVal= K[i] + d[i], interVals=len)}))
  densi2 <- t(sapply(1:nrow(mat), function(i){ Maxent_Density(mat[i,], input, factor, minVal= K[i] + d[i], maxVal= max_limit_int, interVals=len)}))
  d <- d/factor
  K <- K/factor
  integral1 <- NULL
  integral2 <- NULL
  for (j in 1:nrow(densi1)){
    sec  <- seq(d[j], K[j] + d[j],   length = len)
    sec2 <- seq(K[j] + d[j], max_limit_int, length = len)
    
    integral1[j] <- trapz(sec,  (sec - d[j])*densi1[j,]) 
    integral2[j] <-  K[j]*trapz(sec2, densi2[j,]) 
  }
  
  intregal_cap <- (integral1 + integral2)*factor
  #boxplot(integral1 + integral2)
  
  return(intregal_cap)
}

risk_premium <- function(densi, a = 1, type = 'root', factor = 1){
  
  integral <- NULL
  x        <- densi[, 1]
  densi    <- densi[, -1]
  
  if (type == 'exp'){
    #root
    integral <- NULL
    for (j in 1:ncol(densi)){
      #print(i)
      integral[j] <-  -(1/a)*log( trapz(x/factor, densi[,j]*factor*exp(-a*(x/factor) )) )*factor
    }
  }else if(type == 'root'){
    #exp
    for (j in 1:ncol(densi)){
      #print(i)
      integral[j] <-  ((trapz(x/factor, densi[,j]*factor*sqrt(x/factor)))^(2))*factor
    }
  }else{
    cat('Error: type not found')
  }
  
  return(integral)
}

g_distorted_risk_price <- function(densi, CDF, tam, a){
  

  b      <- 1/(1 - exp(-a)) #- 1

  #integral_g     <- sapply(1:nrow(mat), function(i){ round(trapz(x, b*(1-exp(-a*(1 - F_[i,])))), 2)})
  #integral_g     <- t(integral_g)
  
  integral_g <- vector(mode = 'numeric', length = nrow(mat))
  #          g <- matrix(0,  nrow = nrow(mat), ncol = len)
  x <- densi[, 1]
  densi <- densi[, -1]
  CDF <- CDF[, -1]
  #g       <- function(x){b*(1-exp(-a*(1 - x)))}
  g       <- function(x){b*(1-exp(-a*x))}
  g_prima <- function(x){b*(a*exp(-a*x))}
  for (i in 1:ncol(densi)){
  #   if ( all((1 - F_[i,]) == 0)){
  #     g[i, ]        <- 0
  #     integral_g[i] <- 0
  #   }else if( all((1 - F_[i,]) == 1)){
  #     g[i, ] <- 1
  #     integral_g[i] <- 1
  #   }else{
  #     g[i, ] <- b*(1-exp(-a*(1 - F_[i,])))
  #     integral_g[i] <- round(trapz(x, b*(1-exp(-a*(1 - F_[i,])))), 2)
  #   }

     # # integral_g[i] <- trapz(x, g)
     # integral_g[i] <- (F_[i,])*(1-g(1-F_[i])) - trapz(x, 1- g(1-F_[i,]))
     #integral_g[i] <- trapz(x, x*densi1[i,])
     integral_g[i] <- trapz(x, x*(g_prima(1-CDF[,i]))*densi[,i])
     #integral_g[i] <- trapz(x, g(1 - F_[i,]))
   }
  
  
  #saveRDS(integral_g, file =  paste0('input/integral_g_', tam, '.rds'))
  #saveRDS(g, file =  paste0('input/g_', tam, '.rds'))
  #integral_sqrt  <- t(sapply(1:nrow(mat), function(i){trapz(x, sqrt((1 - F_[i,])) )}))
  
  #return(list(g_distorted = integral_g, g_distorted_sqrt = integral_sqrt))
  return(list(g_distorted = integral_g))
}
# 
# ml <- t(rbind(x, F_))
# dits_funct_graph <- function(mat){
#  
#    mat <- as.data.table(cbind(ml[, 1], ml[, 2]))
#    names(mat) <- c('x', 'y')
#   d <- ggplot(mat, aes(x =x, y = y))  + geom_line()
#   for (i in 2:nrow(F_)){
#     mat <- as.data.table(cbind(ml[, 1], ml[, i]))
#     names(mat) <- c('x', 'y')
#     d <- d + geom_line(aes(x =x, y = y))
#   }
#   d <- d + xlab('Sample size') + ylab('stop-loss with a cap') 
#   d <- d + theme(axis.title.y = element_text(size = rel(1.5)))
#   d <- d + theme(axis.text.y = element_text(size = rel(1.8)))
#   d <- d + theme(axis.text.x = element_text(size = rel(1.8)))
#   d <- d + theme(axis.title.x = element_text(size = rel(1.5)))
#   d <- d + theme(legend.text = element_text(size = 15))
#   d <- d + theme(legend.title = element_text(size = 15))
#   d <- d +  theme(
#     panel.background = element_rect(fill="white") ,
#     panel.grid.minor.y = element_line(size=3),
#     panel.grid.major = element_line(colour = "lightgray"),
#     plot.background = element_rect(fill="white")
#   )
#   d+ geom_line()
#   d <- d + xlab('Sample size') + ylab('stop-loss with a cap') 
#   d <- d + theme(axis.title.y = element_text(size = rel(1.5)))
#   d <- d + theme(axis.text.y = element_text(size = rel(1.8)))
#   d <- d + theme(axis.text.x = element_text(size = rel(1.8)))
#   d <- d + theme(axis.title.x = element_text(size = rel(1.5)))
#   d <- d + theme(legend.text = element_text(size = 15))
#   d <- d + theme(legend.title = element_text(size = 15))
#   d <- d +  theme(
#     panel.background = element_rect(fill="white") ,
#     panel.grid.minor.y = element_line(size=3),
#     panel.grid.major = element_line(colour = "lightgray"),
#     plot.background = element_rect(fill="white")
#   )
#   d
# }

