
# input verifications

input_basic_verifications <- function(dt){
  
l = 1
while (l == 1){
  if (length(dt) < 1 | is.null(dt)){
    cat('[Error] dataset is not selected or empty. Please, Choose a valid data set', '\n')
    l = 0
    return()
  }
  if (length(dt)  > 1){
    l = 0
  }
  if (length(dt) < 10 & length(dt)  > 1){
    cat('[Warning] dataset has only a few values. Is this correct?', '\n')
    cat('Number of values: ', length(dt), '\n')
    l = 0
  }
  if (all(is.na(dt))){ #missing values
    cat('[Warning] There are missing values. Erasing missing values:', length(dt[is.na(dt)]), '\n')
    dt <- dt[!is.na(dt)]
    l = 1
  }
  if (!any(sapply(dt, is.numeric))){ # non-numeric values
    cat('[Warning] There are non-numeric values. Erasing non-numeric values:', length(dt[is.na(dt)]), '\n')
    ids <- which(sapply(dt, is.numeric))
    dt  <- dt[, ids]
    l = 1
  }
}
  return(dt)
}


#-----------------------------------------------------------------------------------------------------------------------------


transformation <- function(dt){
  
    num_digits <- nchar(as.character(trunc(max(dt))))
    product    <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits)))))
  # Whe have to scale big numbers to avoid infinite values due to the exponential
  # When you divide you have to take care of the size of the number, small values cannot reach the optimal values
   if (all(dt < 0)) {
     dt <- (-1)*dt
     product <- -1
  }else if (num_digits >= 2 & max(dt) > 20 ){
    if (num_digits < 3){
      dt <- dt/product
    }else{
      num_digits <- num_digits - 1
      product <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits)))))
      dt      <- dt/product
    }
  }else if (max(dt) < 1){
    dt      <- dt*product
  }
    return(list(dt = dt, product = product))
}
#-----------------------------------------------------------------------------------------------------------------------------

apply_transformation_to_other_data <- function(dt, new_data){
  
  num_digits <- nchar(as.character(trunc(max(dt))))
  product    <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits)))))
  # Whe have to scale big numbers to avoid infinite values due to the exponential
  # When you divide you have to take care of the size of the number, small values cannot reach the optimal values
  if (all(dt < 0)) {
    new_data <- (-1)*new_data
    product <- -1
  }else if (num_digits >= 2 & max(dt) > 20 ){
    if (num_digits < 3){
      new_data <- new_data/product
    }else{
      num_digits <- num_digits - 1
      product <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits)))))
      new_data      <- new_data/product
    }
  }else if (max(dt) < 1){
    new_data      <- new_data*product
  }
  return(list(dt = new_data, product = product))
}

#-------------------------------------------------------------------------------------------------------------------------------
#transformation

inverse_transformation <- function(dt_ini, dt, densidad){
  
  num_digits <- nchar(as.character(trunc(max(dt_ini))))
  if (all(dt_ini < 0)) {
    densidad <- densidad %>% mutate(x = (-1)*x, densidad = densidad)
  }else if (num_digits >= 2 & max(dt_ini) > 20){
    if (num_digits < 3){
      product    <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits)))))
      densidad <- densidad %>% mutate(x = x*product, densidad = densidad/product)
    }else{
      product    <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits -1)))))
      densidad <- densidad %>% mutate(x = x*product, densidad = densidad/product)
    }
  } else if (max(dt_ini) < 1){
    product    <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits)))))
    densidad <- densidad %>% mutate(x = x/product, densidad = densidad*product)
  }
  
  return(densidad)
}


F_inverse_transformation <- function(dt_ini, dt, FF){
  
  num_digits <- nchar(as.character(trunc(max(dt_ini))))
  if (all(dt_ini < 0)) {
    FF <- FF %>% mutate(x = (-1)*x)
  }else if (num_digits >= 2 & max(dt_ini) > 20){
    if (num_digits < 3){
      product    <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits)))))
      FF <- FF %>% mutate(x = x*product)
    }else{
      product    <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits -1)))))
      FF <- FF %>% mutate(x = x*product)
    }
  } else if (max(dt_ini) < 1){
    product    <- prod(as.numeric(paste0(1, paste0(rep(0, num_digits)))))
    FF <- FF %>% mutate(x = x/product)
  }
  
  return(FF)
}
