#rLIBARIES
# Required packages list
packages <- c(
  'data.table',
  'caret',
  'plotly',
  'htmlwidgets',
  'htmltools',
  'dplyr',
  'Bolstad',
  'pracma',
  'ggplot2',
  'tidyquant',
  'foreach',
  'highcharter',
  # FitDistributionSimulations
  'fitdistrplus', 
  "actuar",
  'fExtremes',
  'ismev',
  'doParallel'
  
  
)

# sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev


# Install packages not found in local machine
if(any(!(packages %in% installed.packages()))){
  install.packages(packages[!(packages %in% installed.packages())], dependencies = T)
}

# Require all packages
result <- lapply(packages, suppressMessages(require), character.only = TRUE)
