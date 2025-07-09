###################
###combine param estimates and covariates
###################

#load packages that are relevant and a little bit of figure code
source("code/0_libraries.R") 
source("code/depr/functions.R")

#parameters estimated for alpha/beta using nls2
# params<- read.csv("output/params_refined-7-24-23.csv")[,c(2:6)]
# colnames(params)<-c("substudy_num","alpha.nls2","alpha_se.nls2","betanls2","beta_senls2")

params <- read.csv("output/combined_results_2024-09-17.csv")
params <- params[,c(1,3,5)]


#add data that was pre estimate as beta from studies where we didn't fit raw beta estimate
schmitt2007<-read.csv("data/schmitt_holbrook_2007.csv")
schmitt2007$beta_variance<-schmitt2007$se_perm2^2
schmitt2007<-schmitt2007[,c(1,3,7)]
names(schmitt2007)<-names(params)

params <- bind_rows(params,schmitt2007)




wilsonosenberg2002<-read.csv("data/wilson_osenberg_2002.csv")

params <- bind_rows(params,wilsonosenberg2002)



shimaosenberg2003<-read.csv("data/shima_osenberg2003.csv")

params <- bind_rows(params,shimaosenberg2003)


#load covariates
covariates<-read.csv("data/covariates-2024-09-30.csv")

#merge covariates and alpha/beta estimates for just the use = yes studies
covariates<-filter(covariates,use_2024=="yes")
all_dat<-right_join(params,covariates,by="substudy_num")



#pull out study duration by substudy_num and add it to the joined dataset
duration<- read.csv("data/all_studies_looped-2024-09-11.csv")
duration2 <- duration %>%
  dplyr::select(`substudy_num`, `t`)%>%
  group_by(substudy_num) %>%
  dplyr::summarize(duration = mean(t, na.rm=TRUE))

all_dat2<-right_join(duration2,all_dat)


#rename beta to beta_hat
all_dat2 <- rename(all_dat2, beta_hat = beta)


#sample size 2024 - 147
all_dat2%>%
  drop_na(beta_hat)%>%
  tally()



#add in the median density and mean density for each substudy number


# Assuming you have already loaded your data
data <- read.csv("data/all_studies_looped-2024-09-11.csv")
glimpse(data)

# Filter data to retain only matching substudy numbers in all_dat2
matching_data <- data %>%
  filter(substudy_num %in% all_dat2$substudy_num)

# Group by substudy_num and calculate the mean and median of n0_m2
density_stats <- matching_data %>%
  group_by(substudy_num) %>%
  summarise(
    mean_density = mean(n0_m2, na.rm = TRUE),
    median_density = median(n0_m2, na.rm = TRUE)) 

glimpse(density_stats)
# density_stats$study_num <- as.numeric(density_stats$study_num)

#get rid of elactinus empty row 
# density_stats<-density_stats[-72,]



extra_densities = read.csv("data/manual_densities.csv")

#join density stats and extra densities
glimpse(extra_densities)
glimpse(density_stats)



# Merge the data frames
merged_df <- density_stats %>%
  full_join(extra_densities, by = "substudy_num") %>%
  mutate(
    mean_density = coalesce(mean_density.x, mean_density.y)) %>%
  select(substudy_num,mean_density, median_density)  # Select final columns




# Merge the data frames
all_dat2 <- all_dat2 %>%
  left_join(merged_df, by = c("substudy_num"))

write.csv(all_dat2,"output/merged-covariates-10-21-24.csv")




#Number of unique papers
length(unique(all_dat2$study_num))

#number of unique species
length(unique(all_dat2$g_sp))

#variation in duration 
range(all_dat2$duration)




