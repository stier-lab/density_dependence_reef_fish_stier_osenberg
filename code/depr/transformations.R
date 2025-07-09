########################################################################################
##########Data Categorization and Transformations#######################################
########################################################################################

source("code/0_libraries.R")
source("code/4_combine_params_covariates.R")
source("code/1_functions.R")


# #tropical/temperate cutoff
all_dat2$tropical<-ifelse(abs(all_dat2$lat_deci)>30,"Temperate","Tropical")

# #convert from m2 to cm2
all_dat2$betanls2_raw_cm<-all_dat2$beta_hat*1E4
all_dat2$betanlsvar_raw_cm<-all_dat2$beta_variance*(1E4^2)

#converting the beta estimate and variance to asinh units
all_dat2$betanls2_asinh<-asinh(all_dat2$betanls2_raw_cm)
all_dat2$betanlsvar_asinh<-all_dat2$betanlsvar_raw_cm*(1/sqrt((1+(all_dat2$betanlsvar_raw_cm)^2)))


