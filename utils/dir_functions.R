
dir_fig_creation <- function(name.file, method){
  
  # create two levels folder
  fig_dir <- paste0('fig')
  #first level should be created
  if (!dir.exists(fig_dir)){
    cat('[Info]: Creating directory ', fig_dir, '\n')
    dir.create(fig_dir)
  }
  two_step_folder <- c(gsub("\\..*","", name.file), method)
  for (i in two_step_folder){
    fig_dir <- paste0(fig_dir, '/', i )
    if (!dir.exists(fig_dir)){
      cat('[Info]: Creating directory ', fig_dir, '\n')
      dir.create(fig_dir)
    }
  }
  
  return(fig_dir)
}

#--------------------------------------------------------------------------------------

dir_output_creation <- function(name.file, method){
  
  # create two levels folder
  fig_dir <- paste0('output')
  #first level should be created
  if (!dir.exists(fig_dir)){
    cat('[Info]: Creating directory ', fig_dir, '\n')
    dir.create(fig_dir)
  }
  two_step_folder <- c(gsub("\\..*","", name.file), method)
  for (i in two_step_folder){
    fig_dir <- paste0(fig_dir, '/', i )
    if (!dir.exists(fig_dir)){
      cat('[Info]: Creating directory ', fig_dir, '\n')
      dir.create(fig_dir)
    }
  }
  
  return(fig_dir)
}