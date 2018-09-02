

# all the algorithms used in the literature
maxent_types <- function(input){
  
  if(input$method == 'SMEE2'){ # | input$method == 'SMEE2'
    #mp=bootstrap_mu(k,S,seed=10,p1=0.90,p2=1)
    #mp = theorical_interval_mu(ele=1, me=0, sigma=0.25, M=length(S), k=8,rep=10,seed=10) #usado en entropy
    mp = theoretical_interval_mu_new(input$S, k=input$k, M=20, rep=1000, alph=0.05)
    mu_interval=mp #- mu
    write.table(mu_interval, file = "rds/mu_matrix.dat",row.names = FALSE,col.names = FALSE)
  }
  
  name_fun <- 
    switch(input$method,
           SME             = "opt/funSME.R",
           MEM_Poisson     = "opt/funMEMpoisson.R",
           MEM_Exponential = "opt/funMEMexp.R",
           SMEE1           = "opt/fun_mu.R",
           SMEE2           = "opt/fun_mu2.R"
    )
  
  
  name_grad <- 
    switch(input$method,
           SME             = "opt/gradientSME.R",
           MEM_Poisson     = "opt/gradientMEMpoisson.R",
           MEM_Exponential = "opt/gradientMEMexp.R",
           SMEE1           = "opt/grad_mu.R",
           SMEE2           = "opt/grad_mu2.R"
    )
  
  #global <<-1e-3
  #hh <- BB(lambda,alpha,mu,nameFun= name_fun,nameGrad= name_grad,M=1000,tolerance=1e-5,step=1,correction=0)
  
  hh <- BB_modified(input$lambda, input$alpha, input$mu, nameFun = name_fun, nameGrad = name_grad,
                    M = input$iterations,
                    tolerance = input$tol, step = input$step, correction = 0, #input$Correction,
                    input$GLL, sigma1 = input$sigma[1], sigma2 = input$sigma[2],
                    phi_max = input$phi_max1, phi_min = input$phi_min1, N = 1000)
  
  # h <- BBoptim(par = lambda, fn = obj_function_prim, gr =  NULL, method = 3, 
  #              control=list(checkGrad = F, maxit = 10000, #ftol=1.e-10, eps = 1e-07, 
  #                           gtol=1.e-6,  maximize = F)) #gr = gradient_function,  
  
  
  lamb_opt <- hh$lambda_Opt

  return(lamb_opt)
}

generate_x_values <- function(S, num = 300, factor = 2.5){
  sep      <- 15 #set 9 o 10
  d        <- hist(S, breaks = sep, plot = FALSE)
  x        <- seq(0.01, factor*max(d$mids), length = num)
  
  return(x)
} 


density.values <- function(x, input, lamb_opt, globalZ = 1e-3){
  
  alpha <- input$alpha
  if (input$method == "SME"){
    densidad <- SME_density(lamb_opt, alpha, x, globalZ)$densidadS
  }else if (input$method == "MEM_Poisson"){
    densidad <- densityMEMpoisson(lamb_opt, alpha, x, N = 200, ita = 2, meth = 4)
  }else if (input$method == "MEM_Exponential"){
    densidad <- densityMEM_exp(lamb_opt, alpha ,x, N = 190, zeta = 50, mu = NULL, meth = 6)$densityMEM
    #densidad <- densityMEM_exp(lamb_opt, alpha ,x, N = 100, zeta = 10, mu = NULL, meth = 1)$densityMEM
  }else if (input$method == "SMEE1"){
    densidad <- SME_density(lamb_opt, alpha, x, globalZ)$densidadS
  }else if (input$method == "SMEE2"){
    densidad <- SME_density(lamb_opt, alpha, x, globalZ)$densidadS
  }else{
    cat('[INFO] Error')
  }
  return(densidad)
}

