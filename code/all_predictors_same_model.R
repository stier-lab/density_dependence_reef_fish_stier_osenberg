# Load necessary libraries
library(ggplot2)
library(dplyr)
library(metafor)
library(tidyverse)
library(ggcorrplot)
library(corrplot)
library(GGally)
library(plotly)
library(MuMIn)
library(gt)
library(scales)  # For color scaling
library(readxl)
library(metafor)


#####
#Load data and trees
#####

# Load and preprocess data
source("code/4_combine_params_covariates.R")
source("code/0_libraries.R")
source("code/5_transformations.R")


# Load trees
tr <- read.tree("data/Reef_fish_all.tacted.newick.tre")
tr_r <- read.tree("data/actinopt_12k_treePL.tre")
tr_r2 <- read.tree("data/1.newick.tre")
tr_r3 <- fishtree_phylogeny()


# 1. Read in the spreadsheet of max lengths & weights
meta_phy <- read_excel("data/unique_species_studies.xlsx") %>%
  rename(
    max_len_cm = `max_len (cm)`,
    max_wt_g   = `Max_wt (ga)`
  )

# 2. Merge into your main data (all_dat2), by species + study + substudy
all_dat2 <- all_dat2 %>%
  left_join(meta_phy, by = c("g_sp", "study_num", "substudy_num"))

# 3. Compute the new moderator: max length × density
all_dat2 <- all_dat2 %>%
  mutate(
    maxlen_density = max_len_cm * mean_density
  )


#subset list of unique species and their study and substudy number please 

unique_species_studies <- all_dat2 %>%
  select(g_sp, study_num, substudy_num) %>%
  distinct() %>%
  arrange(g_sp, study_num, substudy_num)

# View the result
print(unique_species_studies)

write.csv(unique_species_studies, "data/unique_species_studies.csv", row.names = FALSE)


# Count number of unique species
n_unique_species <- all_dat2 %>%
  distinct(g_sp) %>%
  nrow()

# Print result
cat("Number of unique species:", n_unique_species, "\n")



# 1) Number of studies and sub‐studies
n_studies    <- all_dat2 %>% distinct(study_num)    %>% nrow()
n_substudies <- all_dat2 %>% distinct(substudy_num) %>% nrow()

cat("Number of unique studies:   ", n_studies,    "\n")
cat("Number of unique substudies:", n_substudies, "\n\n")


study_info <- all_dat2 %>%
  select(
    study_num,
    substudy_num,
    Authors,
    Article.Title,
    Source.Title,       # journal name
    Publication.Year
  ) %>%
  distinct() %>%
  arrange(Publication.Year, study_num, substudy_num)

# Quick look
glimpse(study_info)

# Or render as a GT table
study_info %>%
  gt() %>%
  tab_header(
    title    = "Study / Substudy → Paper Details",
    subtitle = "Authors • Article Title • Journal • Year"
  ) %>%
  cols_label(
    study_num       = "Study #",
    substudy_num    = "Substudy #",
    Authors         = "Authors",
    Article.Title   = "Article Title",
    Source.Title    = "Journal",
    Publication.Year = "Year"
  )


library(dplyr)

# 1) Unique study‐substudy rows with publication info
study_substudy_info <- all_dat2 %>%
  select(
    study_num,
    substudy_num,
    Authors,
    Article.Title,
    Source.Title,       # journal
    Publication.Year
  ) %>%
  distinct()

write.csv(
  study_substudy_info,
  file = "output/unique_study_substudy_info.csv",
  row.names = FALSE
)

# 2) Unique study rows (one row per study) with publication info
study_info <- all_dat2 %>%
  select(
    study_num,
    Authors,
    Article.Title,
    Source.Title,
    Publication.Year
  ) %>%
  distinct()

write.csv(
  study_info,
  file = "output/unique_study_info.csv",
  row.names = FALSE
)


########################################################################
##1. Estimate the effect of beta 
######################################################################


# Mixed effects model (Random effects for study and substudy)
beta_ma_all <- rma.mv(yi = betanls2_asinh, V = betanlsvar_asinh, data = all_dat2,
                      random = list(~1|study_num/substudy_num), method = "REML", test = "t")

# Extract coefficients and back-transform them for each model
coef_ma_all <- coef(summary(beta_ma_all))

# Back-transforming the coefficients (using sinh back-transformation and converting cm² to m²)
backtransformed_ma_all <- sinh(coef_ma_all["estimate"]) 
backtransformed_ma_all_lower <- sinh(coef_ma_all["ci.lb"]) 
backtransformed_ma_all_upper <- sinh(coef_ma_all["ci.ub"]) 

# Fixed effect model: extract total tau² (heterogeneity)
# tau2_fixed <- summary(beta_ma_all)$tau2  # Total tau² from the fixed effects model
# tau_sd <- sqrt(tau2_fixed)
# tau_bt<- sinh(tau_sd)


#these are on transformed scale 
# Nested random effects model: extract study and substudy-level variances
tau2_study <- summary(beta_ma_all)$sigma2[1]   # Variance for study level
tau2_substudy_nested <- summary(beta_ma_all)$sigma2[2]  # Variance for substudy level


# The variance components extracted from summary(beta_ma_all)$sigma2 are in the arcsinh-transformed scale.
# To approximate these variances in the original scale, we use the Delta method.

# The Delta method states that for a transformed variable Y = f(X), 
# the variance of Y can be approximated as: Var(Y) ≈ [df/dX]^2 * Var(X),
# where df/dX is the derivative of the transformation function.

# Since we used the arcsinh transformation, the inverse function is sinh(X).
# The derivative of sinh(X) is cosh(X), so we approximate the variance in the original scale as:

# Var(original) ≈ sinh(sqrt(Var(transformed)))^2

# Apply this transformation to back-transform the variance components:

backtransformed_tau2_study <- sinh(sqrt(tau2_study))^2
backtransformed_tau2_substudy <- sinh(sqrt(tau2_substudy_nested))^2

# Print the results
cat("Back-transformed study-level variance:", backtransformed_tau2_study, "\n")
cat("Back-transformed substudy-level variance:", backtransformed_tau2_substudy, "\n")





######################################################################
#2. Pull out species estimates and map on phylogeny
######################################################################


# Process the all_clean data
all <- all_dat2 %>%
  select(g_sp, betanls2_raw_cm,betanls2_asinh, betanlsvar_asinh,
         study_num, substudy_num, family, predators,use_2024,expt_obs,
         size_start, expt_obs_pairs,duration, mean_density,family,paired_pred,paired_substudy_num,max_len_cm,max_wt_g)
#%>%
  # filter(g_sp != "Sebastes_spp_NA") %>%
 # filter(predators=="present")

# Get unique species from all_dat and all
all_dat_species <- unique(all_dat2$g_sp)
all_species <- unique(all$g_sp)

# Compare the lengths
length_all_dat <- length(all_dat_species)
length_all <- length(all_species)

# Find species in all_dat but not in all
missing_in_all <- setdiff(all_dat_species, all_species)

# Find species in all but not in all_dat
missing_in_all_dat <- setdiff(all_species, all_dat_species)

# Print results
cat("Number of unique species in all_dat:", length_all_dat, "\n")
cat("Number of unique species in all:", length_all, "\n")
cat("Species in all_dat but not in all:\n", missing_in_all, "\n")
cat("Species in all but not in all_dat:\n", missing_in_all_dat, "\n")

#Epinephelus_malabaricus is from Chua and Teng - an aquaculture study
#Fundulus_heteroclitus is from Kneib 

# Replace  "Sebastes_spp" with "Sebastes_mystinus" in the g_sp column of all
#assuming sebastes is far enough away that the KGB group all can be one one tip 

#change some species names
all$g_sp <- gsub("Sebastes_spp_NA", "Sebastes_mystinus", all$g_sp)
all$g_sp <- gsub("Fundulus_heteroclitus", "Fundulus_heteroclitus_heteroclitus", all$g_sp)
all$g_sp <- gsub("Diplodus_sargus", "Diplodus_sargus_sargus", all$g_sp)
all$g_sp <- gsub("Elactinus_sp\\.", "Elacatinus_evelynae", all$g_sp)

sort(unique(all$g_sp))

#note that fundulus and epineph are both predators = absent

# Raw means for each species
beta_raw <- all %>%
  group_by(g_sp) %>%
  summarise(mean = mean(betanls2_asinh))

# Filter for non-missing values
all_clean <- all %>%
  filter(!is.na(betanls2_asinh) & !is.na(betanlsvar_asinh))

# Filter the all_clean data to include only species in species_counts
filtered_data <- all 



##estimates of beta for phylogeny##

#look at species effects, ignoring phylo covariance
m_sp <- rma.mv(yi = betanls2_asinh, V = betanlsvar_asinh, data = filtered_data, 
                  mods = ~ g_sp,
                  random = list(~ 1 | study_num / substudy_num), 
                  method = "REML", test = "t")

summary(m_sp)

# Run the random effects model without intercept for species-specific estimates 
#no phylo covariance here because we're trying 

m_sp_ni <- rma.mv(yi = betanls2_asinh, V = betanlsvar_asinh, data = filtered_data, 
                  mods = ~ g_sp - 1,  # No intercept
                  random = list(~ 1 | study_num / substudy_num), 
                  method = "REML", test = "t")


m_sp_ni_coef<-summary(m_sp_ni)$b[,1]

###---Pull out coefficient estimates and compare to raw beta means---####

model_species <- rownames(coef(summary(m_sp_ni)))
print("Model species:")
print(model_species)

# Remove the prefix from model_species
clean_model_species <- gsub("g_sp", "", model_species)

# Print the cleaned model species
print("Cleaned model species:")
print(clean_model_species)

# Identify species present in the model and in the data
data_species <- sort(unique(all_dat2$g_sp))

# Check for missing species
missing_species <- setdiff(data_species, clean_model_species)
print("Missing species:")
print(missing_species)

#note that fundulus and epineph are both predators = absent
#sebastes_spp is multi speceies so excluded

# Create the data frame with aligned species data
beta_sp_df <- data.frame(
  # "betaraw" = sinh(beta_raw_filtered$mean),
  "beta_metafor_fe" = sinh(m_sp_ni_coef),
  "g_sp" = clean_model_species
)

print(beta_sp_df)
# 
# #for this phylogeny  # Replace "Elactinus_sp." with "Elacatinus_evelynae" in the g_sp column of beta_sp_df
# beta_sp_df$g_sp <- gsub("Elactinus_sp\\.", "Elacatinus_evelynae", beta_sp_df$g_sp)
# 
# #rplace Fundulus_heteroclitus with Fundulus_heteroclitus_heteroclitus
# beta_sp_df$g_sp <- gsub("Fundulus_heteroclitus", "Fundulus_heteroclitus_heteroclitus", beta_sp_df$g_sp)
# 
# #replace diplodus sargus with Diplodus_sargus_sargus 
# beta_sp_df$g_sp <- gsub("Diplodus_sargus", "Diplodus_sargus_sargus", beta_sp_df$g_sp)
# Verify the replacement

##---Compare the meta analysis and tree dataset species list ---##

# unique species list from our data
our_sp_unique <- beta_sp_df$g_sp

# Which species are in the phylogeny?
tree_sp <- TipLabels(tr_r2)

d1 <- data.frame("genus_species" = tree_sp)
d2 <- data.frame("genus_species" = our_sp_unique)

# Check species present in both datasets
matching_species <- intersect(d1$genus_species, d2$genus_species) 
missing_in_tree <- setdiff(d2$genus_species, d1$genus_species)
missing_in_data <- setdiff(d1$genus_species, d2$genus_species)

print("Species in data but not in phylogeny:")
print(missing_in_tree)



#Diplodus_sargus not in tree

# print("Species in phylogeny but not in data:")
# print(missing_in_data)

#three species aren't present, diplodus sargus, elactinus and the sebastes_spp 

# How many species do we have?
num_species <- length(our_sp_unique) 
print(paste("Number of unique species in our data:", num_species))

# List of matching species in the phylogeny
matchy <- intersect(tree_sp, our_sp_unique)
print("Matching species in phylogeny:")
print(matchy)

# Look for misspellings or case differences
species_mismatch <- stringdist_inner_join(d1, d2, ignore_case = FALSE, distance_col = "distance")
print("Species with possible misspellings or case differences:")
print(species_mismatch)

##---Adjust names for two species with intraspecies identity---## 

# names(beta2)<-str_replace(rownames(beta2), "Fundulus_heteroclitus", "Fundulus_heteroclitus_heteroclitus")
# names(beta2)<-str_replace(names(beta2),"Diplodus_sargus", "Diplodus_sargus_sargus")

# betadf2$g_sp<-str_replace(betadf2$g_sp, "Fundulus_heteroclitus", "Fundulus_heteroclitus_heteroclitus")
# beta_sp_df$g_sp<-str_replace(beta_sp_df$g_sp, "Diplodus_sargus", "Diplodus_sargus_sargus")


#---Prune tree to just species we have data from rabowski paper---#

#prune tree to just species we have data from rabowski paper 

beta_sp_df<-
  beta_sp_df%>%
  filter(g_sp != "Diplodus_sargus")

beta2 <- beta_sp_df$beta_metafor_fe
names(beta2) <- beta_sp_df$g_sp

# Ensure names consistency and print mismatches
matched_species <- intersect(names(beta2), tree_sp)
if (length(matched_species) != length(beta2)) {
  warning("Some species in beta2 do not match species in the tree.")
  print(setdiff(names(beta2), tree_sp))
}

# Prune the tree to include only species we have data for
tr4 <- drop.tip(tr_r2, setdiff(tr_r2$tip.label, matched_species))

TipLabels(tr4)


# Ensure names match between beta2 and tree_sp
matched_species <- intersect(names(beta2), tree_sp)
if (length(matched_species) != length(beta2)) {
  warning("Some species in beta2 do not match species in the tree.")
  print(setdiff(names(beta2), tree_sp))
}


betavec<-c(beta)
names(betavec)<-rownames(beta)
# tr5<-prune.missing(betavec,tr_r)

head(beta_sp_df)

#log scale beta plot with tips as meta analysis estimates , swaw btw betaraw and beta_metafor
beta_sp_df2<-beta_sp_df%>%
  select(g_sp,"mean"=beta_metafor_fe)

ggtree(tr4,layout="circular")%<+%
  beta_sp_df2+
  geom_tippoint(aes(colour=mean),size=3)+
  scale_colour_gradientn(colours = pal,name="Beta",trans='asinh',
                         breaks=c(-10,0,10,100,1000,10000))+
  geom_tiplab(size=2,offset=2)


range(beta_sp_df$beta_metafor_fe)

# Plot the circular phylogeny
p <- ggtree(tr4, layout = "circular") %<+% beta_sp_df2 +
  
  # First layer for black outline of points
  geom_tippoint(aes(fill = mean), shape = 21, size = 5, color = "black", stroke = 1) +
  
  # Gradient fill for continuous mean values (white to red)
  scale_fill_gradient(low = "white", high = "red", 
                      name = "Beta", trans = 'asinh',
                      limits = c(0, 10000),  # Ensure the scale covers 0 to 10,000
                      breaks = c(0, 10, 100, 1000, 10000),
                      labels = c("0", "10", "100", "1,000", "10,000"),
                      guide = guide_colorbar(frame.colour = "black", 
                                             frame.linewidth = .5,
                                             ticks.colour = "black", 
                                             ticks.linewidth = .5)) +
  # Add tip labels
  geom_tiplab(size = 2, offset = 2) +
  
  # Tree theme
  theme_tree2()+
  theme(
    axis.line = element_blank(),        # Remove axis lines
    panel.grid = element_blank(),        # Remove any panel grid lines
    plot.margin = margin(10, 10, 10, 10) # Optional: Adjust plot margins
  )

# Display the plot
print(p)



###-- estimate whether beta is nonrandom bloomberg's k pagel's lambda-- ### 

#create vectof mean beta by species
betavec2<-beta_sp_df2$mean
names(betavec2)<-beta_sp_df2$g_sp

#pagel's lambda 
name.check(tr4,betavec2)

#pagel's lambda with phytools 

#If our estimated lambda = 0, then the traits are inferred to have
#no phylogenetic signal (i.e., traits evolve independent of phylogeny)

# pagels lambda
plot.phylosig(phylosig(tr4,betavec2,method="lambda",test=TRUE)) 
# Phylogenetic signal lambda : 0.637898 
# logL(lambda) : -318.857 
# LR(lambda=0) : 12.6645 
# P-value (based on LR test) : 0.000372659 

#bloomberg's K
phylosig(tr4,betavec2,method="K",test=TRUE) #K = 0.540141 , P-value (based on 1000 randomizations) : 0.1 













##############################################################################################
### Predation effects on density dependence ###
##############################################################################################

#fraction of studies predator present 

# Count the number of studies where predators are present vs. absent
predator_counts <- all_dat2 %>%
  group_by(predators) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))  # Compute fraction

# Print results
print(predator_counts)

# Alternatively, print formatted output
cat("Fraction of studies with predators present:", round(predator_counts$proportion[predator_counts$predators == "present"], 3), "\n")
cat("Fraction of studies with predators absent:", round(predator_counts$proportion[predator_counts$predators == "absent"], 3), "\n")

#make a variance covariance matrix for the phylogeny based on the tree 
phylo_vcv <- vcv(tr4, corr = TRUE)  # Use 'corr = TRUE' to standardize to correlations


# Extract species from dataset
species_in_data <- sort(unique(all_dat2$g_sp))

# Extract species in the phylogenetic matrix
species_in_phylo <- sort(rownames(phylo_vcv))

# Identify species in data but missing in phylogeny
missing_species <- setdiff(species_in_data, species_in_phylo)

# Print missing species
cat("Species in dataset but not in phylogenetic matrix:\n")
print(missing_species)


#change some species names - only if above isn't run 
# all$g_sp <- gsub("Sebastes_spp_NA", "Sebastes_mystinus", all$g_sp)
# all$g_sp <- gsub("Fundulus_heteroclitus", "Fundulus_heteroclitus_heteroclitus", all$g_sp)
# all$g_sp <- gsub("Diplodus_sargus", "Diplodus_sargus_sargus", all$g_sp)
# all$g_sp <- gsub("Elactinus_sp\\.", "Elacatinus_evelynae", all$g_sp)

# # Extract species from dataset
# species_in_data <- sort(unique(all$g_sp))
# 
# # Identify species in data but missing in phylogeny
# missing_species <- setdiff(species_in_data, species_in_phylo)
# 
# # Print missing species
# cat("Species in dataset but not in phylogenetic matrix:\n")
# print(missing_species)
# 
# 
# # Filter dataset to only include species present in phylo_vcv
# all_filtered <- all[all$g_sp %in% species_in_phylo, ]
# 
# # Extract species names from the dataset and phylogenetic matrix
# species_in_data <- unique(all_dat2_filtered$g_sp)  # Only species present in filtered dataset
# species_in_phylo <- rownames(phylo_vcv)  # All species in phylogenetic matrix
# 
# # Find species in the phylogenetic matrix but NOT in the dataset
# extra_species <- setdiff(species_in_phylo, species_in_data)
# 
# # Print extra species
# cat("Species in phylogenetic matrix but not in dataset:\n")
# print(extra_species)

# Re-run the model with the cleaned dataset
m_pred <- rma.mv(yi = betanls2_asinh, 
                 V = betanlsvar_asinh, 
                 data = all,
                 mods = ~predators,
                 random = list(~1 | study_num/substudy_num,
                               ~1 | g_sp),
                 R = list(g_sp = phylo_vcv),
                 method = "REML",
                 test = "t")

summary(m_pred)



# Extract coefficients from the model summary
model_coeffs <- coef(summary(m_pred))

# Get effect size for predators absent (intercept)
effect_no_predators <- model_coeffs["intrcpt", "estimate"]

# Get effect size for predators present (intercept + predator effect)
effect_with_predators <- effect_no_predators + model_coeffs["predatorspresent", "estimate"]

# Compute percentage change
percent_change <- ((effect_with_predators - effect_no_predators) / effect_no_predators) * 100

# Print results
cat("Effect size without predators:", round(effect_no_predators, 4), "\n")
cat("Effect size with predators:", round(effect_with_predators, 4), "\n")
cat("Percentage change in effect size due to predator presence:", round(percent_change, 2), "%\n")


# Extract coefficient estimates
beta_predators <- m_pred$b["predatorspresent", 1]  # Effect of predators
beta_no_predators <- m_pred$b["intrcpt", 1]  # Effect when predators are absent

# Extract confidence intervals (correct indexing)
ci_no_predators <- c(m_pred$ci.lb[1], m_pred$ci.ub[1])  # First row corresponds to "intrcpt"
ci_predators <- c(m_pred$ci.lb[2], m_pred$ci.ub[2])  # Second row corresponds to "predatorspresent"

# Back-transform using sinh() function
beta_no_predators_orig <- sinh(beta_no_predators)
beta_predators_orig <- sinh(beta_no_predators + beta_predators)  # Apply moderator effect

# Back-transform confidence intervals
ci_no_predators_orig <- sinh(ci_no_predators)  # Apply sinh to the CIs
ci_predators_orig <- sinh(ci_no_predators + ci_predators)  # Apply sinh transformation

# Print results
cat("Effect size without predators: β =", round(beta_no_predators_orig, 3), 
    ", CI:", round(ci_no_predators_orig[1], 3), "-", round(ci_no_predators_orig[2], 3), "\n")

cat("Effect size with predators: β =", round(beta_predators_orig, 3), 
    ", CI:", round(ci_predators_orig[1], 3), "-", round(ci_predators_orig[2], 3), "\n")


# Create a data frame with key results
results_table <- data.frame(
  Metric = c("Studies with predators present",
             "Studies with predators absent",
             "Fraction with predators present",
             "Fraction with predators absent",
             "Effect size without predators (β)",
             "Effect size with predators (β)",
             "Percentage change in effect size"),
  
  Estimate = c(predator_counts$count[predator_counts$predators == "present"],
               predator_counts$count[predator_counts$predators == "absent"],
               round(predator_counts$proportion[predator_counts$predators == "present"], 3),
               round(predator_counts$proportion[predator_counts$predators == "absent"], 3),
               round(beta_no_predators_orig, 3),
               round(beta_predators_orig, 3),
               paste0(round(percent_change, 2), "%")), # Format percentage nicely
  
  CI = c("",
         "",
         "",
         "",
         paste0("[", round(ci_no_predators_orig[1], 3), ", ", round(ci_no_predators_orig[2], 3), "]"),
         paste0("[", round(ci_predators_orig[1], 3), ", ", round(ci_predators_orig[2], 3), "]"),
         "")
)

# Generate a high-quality gt table
results_table %>%
  gt() %>%
  tab_header(
    title = md("**Summary of Key Results**"),
    subtitle = "Effect of Predator Presence on Study Outcomes"
  ) %>%
  fmt_number(columns = "Estimate", decimals = 3) %>%
  cols_label(
    Metric = "Metric",
    Estimate = "Estimate",
    CI = "95% Confidence Interval"  # Corrected column name reference
  ) %>%
  tab_options(
    table.font.size = px(14),
    table.border.top.style = "solid",
    table.border.top.color = "black",
    table.border.top.width = px(2),
    table.border.bottom.style = "solid",
    table.border.bottom.color = "black",
    table.border.bottom.width = px(2),
    column_labels.font.weight = "bold",
    heading.align = "center"
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(rows = 5:6) # Bold effect size rows
  )

# ──────────────────────────────────────────────────────────────────────
# Compare Predator Effects: All Studies vs Paired Studies (with low/high % change)
# ──────────────────────────────────────────────────────────────────────

library(metafor)
library(dplyr)
library(tibble)
library(gt)

# Helper to summarize predator effect
summarize_pred_effect <- function(model) {
  cfs      <- coef(summary(model))
  int_log  <- cfs["intrcpt",          "estimate"]
  mod_log  <- cfs["predatorspresent", "estimate"]
  
  # back‐transform
  beta_no    <- sinh(int_log)
  beta_yes   <- sinh(int_log + mod_log)
  pct_change <- 100 * (beta_yes - beta_no) / abs(beta_no)
  
  # CIs on log scale
  ci_int <- c(model$ci.lb[1], model$ci.ub[1])
  ci_mod <- c(model$ci.lb[2], model$ci.ub[2])
  
  # back‐transform CIs
  ci_no  <- sinh(ci_int)
  ci_yes <- sinh(int_log + ci_mod)
  
  # approximate percent‐change bounds
  pct_low  <- 100 * (ci_yes[1] - ci_no[2]) / abs(ci_no[2])
  pct_high <- 100 * (ci_yes[2] - ci_no[1]) / abs(ci_no[1])
  
  tibble(
    beta_no      = beta_no,
    ci_no_low    = ci_no[1],
    ci_no_high   = ci_no[2],
    beta_yes     = beta_yes,
    ci_yes_low   = ci_yes[1],
    ci_yes_high  = ci_yes[2],
    pct_mean     = pct_change,
    pct_low      = pct_low,
    pct_high     = pct_high
  )
}

# Fit on all studies
m_all <- rma.mv(
  yi    = betanls2_asinh,
  V     = betanlsvar_asinh,
  data  = all,
  mods  = ~ predators,
  random= list(~1 | study_num/substudy_num, ~1 | g_sp),
  R     = list(g_sp = phylo_vcv),
  method= "REML", test = "t"
)

# Fit on paired studies only
paired_df <- all %>% filter(paired_pred == "paired")
m_paired <- rma.mv(
  yi    = betanls2_asinh,
  V     = betanlsvar_asinh,
  data  = paired_df,
  mods  = ~ predators,
  random= list(~1 | study_num/substudy_num, ~1 | g_sp),
  method= "REML", test = "t"
)

# Summarize both
sum_all    <- summarize_pred_effect(m_all)    %>% mutate(set = "All studies")
sum_paired <- summarize_pred_effect(m_paired) %>% mutate(set = "Paired only")

# Bind & format
out_df <- bind_rows(sum_all, sum_paired) %>%
  transmute(
    `Study set`            = set,
    `β no predators`       = sprintf("%.2f [%.2f, %.2f]",
                                     beta_no, ci_no_low, ci_no_high),
    `β with predators`     = sprintf("%.2f [%.2f, %.2f]",
                                     beta_yes, ci_yes_low, ci_yes_high),
    `% increase (mean)`    = sprintf("%.1f%%", pct_mean),
    `% inc (low bound)`    = sprintf("%.1f%%", pct_low),
    `% inc (high bound)`   = sprintf("%.1f%%", pct_high)
  )

# Display with gt
out_df %>%
  gt() %>%
  tab_header(
    title    = md("**Predator Presence Effects on Density-Dependence**"),
    subtitle = "All studies vs. Paired-only subset"
  ) %>%
  cols_label(
    `Study set`          = "Study set",
    `β no predators`     = md("**β** (no predators)"),
    `β with predators`   = md("**β** (with predators)"),
    `% increase (mean)`  = md("**% inc** (mean)"),
    `% inc (low bound)`  = md("**% inc** (low)"),
    `% inc (high bound)` = md("**% inc** (high)")
  ) %>%
  tab_options(
    table.font.size = px(14),
    heading.align   = "center"
  )

# Create a dataset for predictions
new_pred_data_pred <- data.frame(predators = factor(c("absent", "present"), levels = c("absent", "present")))

# Generate model matrix for `predators`
model_matrix_pred <- model.matrix(~ predators, data = new_pred_data_pred)

# Extract the dummy variable for "present" (since "absent" is baseline)
predators_dummy <- model_matrix_pred[, "predatorspresent", drop = FALSE]

# Run predictions using `m_pred`
predicted_values_pred <- predict(m_pred, newmods = as.matrix(predators_dummy))

# Store predictions
new_pred_data_pred$predicted_sinhbeta <- predicted_values_pred$pred
new_pred_data_pred$predicted_beta <- sinh(new_pred_data_pred$predicted_sinhbeta)

# Print final prediction dataset
print(new_pred_data_pred)

# Define colors
pred_absent_color <- "#1E90FF"  # Dodger Blue for no predators
pred_present_color <- "#FF4500"  # Red-Orange for predator presence

# Create the plot
ggplot(all, aes(x = predators, y = betanls2_raw_cm, fill = predators)) +
  
  # Boxplots to summarize data distribution
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5, color = "black") +
  
  # Raw data points with slight jitter for visibility
  geom_jitter(aes(color = predators), width = 0.2, alpha = 0.6, size = 3) +
  
  # Predicted points (black dots)
  geom_point(data = new_pred_data_pred, aes(x = predators, y = predicted_beta), 
             shape = 21, size = 5, fill = "black", color = "black", stroke = 1.5) +
  
  # Dashed horizontal reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  
  # Y-axis transformed for better spread of values
  scale_y_continuous(
    trans = 'asinh',
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000, 10000)
  ) +
  
  # Custom colors for predator presence vs. absence
  scale_fill_manual(values = c("absent" = pred_absent_color, "present" = pred_present_color)) +
  scale_color_manual(values = c("absent" = pred_absent_color, "present" = pred_present_color)) +
  
  # Axis labels
  labs(
    x = "Predator Presence",
    y = "Strength of Density-Dependent Mortality (β)"
  ) +
  
  # Clean theme
  theme_classic(base_size = 16) +  
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "gray90"),  
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save the plot
ggsave(
  filename = "figures/predators_vs_beta_boxplot.png",
  width = 5,    
  height = 5,   
  units = "in",
  dpi = 300       
)

# Print final predictions
print(new_pred_data_pred)



#include paired data 


# Define colors for the updated aesthetic
steel_blue <- "#4682B4"  # Used for beta points
dark_blue <- "#1E3A5F"  # Used for outlines of points
paired_color1 <- "#D55E00"  # Complementary orange for paired studies
paired_color2 <- "#009E73"  # Complementary green for paired studies

# Ensure predator labels are capitalized for consistency
all_pred <- all %>%
  mutate(predators = toupper(predators))  # Convert 'absent'/'present' to 'ABSENT'/'PRESENT'

# Separate paired and unpaired studies
paired_studies <- all_pred %>%
  filter(expt_obs == "Exp", paired_pred == "paired")

unpaired_studies <- all_pred %>%
  filter(expt_obs == "Exp", paired_pred != "paired")

# Ensure `paired_substudy_num` exists, assign if missing
if (!"paired_substudy_num" %in% colnames(paired_studies)) {
  paired_studies <- paired_studies %>%
    mutate(paired_substudy_num = as.factor(substudy_num))  
}

# Convert to factors for proper grouping
paired_studies$paired_substudy_num <- as.factor(paired_studies$paired_substudy_num)
unpaired_studies$paired_substudy_num <- as.factor(unpaired_studies$substudy_num)

# Generate distinct colors for paired studies using the updated complementary colors
num_colors <- length(unique(paired_studies$paired_substudy_num))
paired_colors <- setNames(colorRampPalette(c(paired_color1, paired_color2))(num_colors), 
                          unique(paired_studies$paired_substudy_num))

# Create dataset for meta-analysis predictions
new_pred_data_pred <- expand.grid(predators = c("ABSENT", "PRESENT"))

# Ensure predators is a factor
new_pred_data_pred$predators <- factor(new_pred_data_pred$predators, levels = c("ABSENT", "PRESENT"))

# Generate model matrix for `predators`
model_matrix_pred <- model.matrix(~ predators, data = new_pred_data_pred)

# Extract the dummy variable for "PRESENT" (since "ABSENT" is baseline)
predators_dummy <- model_matrix_pred[, "predatorsPRESENT", drop = FALSE]

# Run predictions using `m_pred`
predicted_values_pred <- predict(m_pred, newmods = as.matrix(predators_dummy))

# Store predictions with confidence intervals
new_pred_data_pred <- new_pred_data_pred %>%
  mutate(
    predicted_sinhbeta = predicted_values_pred$pred,
    beta_hat = sinh(predicted_sinhbeta),
    lower = sinh(predicted_values_pred$ci.lb),
    upper = sinh(predicted_values_pred$ci.ub)
  )



# Generate the plot with the improved aesthetic
g <- ggplot() +
  
  # Unpaired studies (jittered points)
  geom_jitter(data = unpaired_studies, 
              aes(x = predators, y = sinh(betanls2_asinh), fill = predators), 
              alpha = 1, height = 0, width = 0.1, 
              shape = 21, color = dark_blue, 
              size = 2.5, stroke = 1) +
  
  # Connect paired studies with geom_path()
  geom_path(data = paired_studies, 
            aes(x = predators, y = sinh(betanls2_asinh), 
                group = paired_substudy_num, 
                color = paired_substudy_num), 
            linewidth = 1, alpha = 0.9) +
  
  # Overlay paired studies with points
  geom_point(data = paired_studies, 
             aes(x = predators, y = sinh(betanls2_asinh), 
                 color = paired_substudy_num, 
                 fill = paired_substudy_num), 
             shape = 21, size = 3, stroke = 1) + 
  
  # Error bars for meta-analysis results
  geom_errorbar(data = new_pred_data_pred, 
                aes(x = predators, y = beta_hat, ymin = lower, ymax = upper), 
                width = 0.05, size = 1, color = dark_blue) +
  
  # Summary points for meta-analysis results
  geom_point(data = new_pred_data_pred, aes(x = predators, y = beta_hat), 
             shape = 18, size = 5, color = dark_blue) +
  
  # Custom color mapping:
  scale_color_manual(values = paired_colors) +
  
  # Fill colors: 
  #   - Unpaired studies → "white" for ABSENT, "steel_blue" for PRESENT
  #   - Paired studies → their respective group color
  scale_fill_manual(values = c("ABSENT" = "white", "PRESENT" = steel_blue, paired_colors)) +
  
  # Styling
  theme_classic(base_size = 12) +
  
  # Y-axis transformation
  scale_y_continuous(trans = 'asinh', 
                     breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)) +
  
  # X-axis labels
  scale_x_discrete(labels = c("Absent", "Present")) + 
  
  # Horizontal reference line at y=0
  geom_hline(yintercept = 0, color = "darkgray", linetype = "dashed") +
  
  # Labels
  labs(
    x = "Predator Treatment",
    y = "Strength of Density-Dependent Mortality (β)"
  ) +
  
  # Text adjustments for a publication-ready plot
  theme(
    legend.position = "none",  
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks.length = unit(0.15, "cm"),
    plot.margin = margin(5, 5, 5, 5)  
  )

# Print the final plot
print(g)

# Save the plot with high resolution
ggsave(
  filename = "figures/predators_vs_beta_paired.png",
  width = 6,  
  height = 6,  
  units = "in",
  dpi = 300      
)







##################################################################################################
#Autocorrelation among predictor variables  
##################################################################################################


all$logmeandensity<-log(all$mean_density)


# Select relevant predictor variables
predictors <- all %>%
  select(expt_obs,size_start, duration, logmeandensity,max_len_cm,max_wt_g)

# Create a multipanel scatterplot matrix
ggpairs(predictors, 
        lower = list(continuous = wrap("smooth", color = "blue")),  # Add smoothed regression lines
        diag = list(continuous = "densityDiag"),  # Use density plots on the diagonal
        upper = list(continuous = wrap("cor", size = 5))) +  # Display correlation coefficients
  theme_minimal()


#look at how much density varies with study type
# Define log-scaled tick marks for density
log_ticks <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)

# Define sinh-scaled tick marks for beta
sinh_ticks <- sinh(c(-10, -5, -2, -1, 0, 1, 2, 5, 10))

ggplot(all, aes(x = mean_density, y = betanls2_raw_cm, shape = expt_obs, color = expt_obs)) +
  geom_point(alpha = 0.5, size = 2) +  # Scatter points with different symbols
   geom_abline(intercept = 0, slope = 0, linetype = "dashed", color = "black") +  # Zero reference line
  scale_x_continuous(
    breaks = log_ticks,
    labels = log_ticks,
    trans = "log10"
  ) +  # Log-scaled x-axis for density
  scale_y_continuous(
    trans = "asinh", 
    breaks = sinh_ticks
  ) +  # Sinh-transformed y-axis for beta
  scale_shape_manual(values = c(16, 17)) +  # Assigns different point shapes for expt/obs
  scale_color_manual(values = c("blue", "red")) +  # Custom colors for expt/obs
  labs(
    x = "Mean Density (n0_m2) by Study (Log Scale)",
    y = "Strength of Density Dependent Mortality (Beta)",
    title = "Beta vs. Density by Study Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "top"
  )












##################################################################################################
#Full model and stepAIC to identify the best fit model 
##################################################################################################

#we look at main effects 

#pull vcov matrix from the tree, this is going to missing or mismatch af a few species names 

phylo_vcv <- vcv(tr4, corr = TRUE)  # Use 'corr = TRUE' to standardize to correlations

all <- all %>%
  filter(predators=="present")

all$logmeandensity<-log(all$mean_density)

# Get unique species in the dataset
species_in_data <- sort(unique(all$g_sp))

# Get species in the phylogenetic matrix
species_in_phylo <- rownames(phylo_vcv)

# Find species present in phylo_vcv but missing in all
setdiff(species_in_phylo, species_in_data)

#need to drop both of these species from the phylogenetic matrix because they're pred absent only
# Find species that are in phylo_vcv but NOT in all
species_to_remove <- c("Epinephelus_malabaricus", "Fundulus_heteroclitus_heteroclitus")

# Keep only species that are in both datasets
species_to_keep <- setdiff(rownames(phylo_vcv), species_to_remove)

# Subset the phylogenetic matrix
phylo_vcv <- phylo_vcv[species_to_keep, species_to_keep]



m_all<-rma.mv(yi=betanls2_asinh, V=betanlsvar_asinh, data=all,
              mods=~expt_obs*logmeandensity*duration*size_start*max_len_cm,
              random = list(~1 | study_num/substudy_num, ~1 | g_sp),
              R = list(g_sp = phylo_vcv),
              method="REML",test="t")

summary(m_all)





# Define the global model predictors
predictors <- c("expt_obs", "size_start", "duration", "logmeandensity","max_len_cm")

# Generate all possible combinations of predictors
all_combinations <- unlist(lapply(1:length(predictors), 
                                  function(x) combn(predictors, x, simplify = FALSE)), recursive = FALSE)

# Placeholder for results
model_results <- data.frame(
  Model = character(),
  Predictors = character(),
  logLik = numeric(),
  AIC = numeric(),
  AICc = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each combination of predictors and fit models
for (comb in all_combinations) {
  formula <- as.formula(paste("~", paste(comb, collapse = " + ")))
  
  # Fit the model with the current combination of predictors
  result <- tryCatch({
    model <- rma.mv(yi = betanls2_asinh, V = betanlsvar_asinh, 
                    mods = formula,
                    random = list(~1 | study_num/substudy_num, ~1 | g_sp),
                    R = list(g_sp = phylo_vcv), 
                    data = all, 
                    method = "REML")
    
    # Extract model summary
    model_summary <- summary(model)
    
    # Extract Log-Likelihood (ll), AIC, and AICc from fit.stats
    logLik_value <- model_summary$fit.stats["ll", "REML"]  # Correct extraction of logLik
    aic_value <- model_summary$fit.stats["AIC", "REML"]    # Extract AIC
    aicc_value <- model_summary$fit.stats["AICc", "REML"]  # Extract AICc
    
    # Return a properly structured data frame
    data.frame(
      Model = paste(comb, collapse = " + "),
      Predictors = paste(comb, collapse = ", "),
      logLik = logLik_value,
      AIC = aic_value,
      AICc = aicc_value
    )
  }, error = function(e) {
    # Return a single-row data frame with NA values in case of failure
    data.frame(
      Model = paste(comb, collapse = " + "),
      Predictors = paste(comb, collapse = ", "),
      logLik = NA,
      AIC = NA,
      AICc = NA
    )
  })
  
  # Append results to the dataframe
  model_results <- rbind(model_results, result)
}

# View final model results
print(model_results)



# Sort model results by AICc in ascending order
model_results_sorted <- model_results %>%
  arrange(AICc)

# Create a nicely formatted table with color formatting
model_results_table <- model_results_sorted %>%
  gt() %>%
  tab_header(
    title = "Model Selection Results",
    subtitle = "Sorted by AICc (Lower is Better)"
  ) %>%
  cols_label(
    Model = "Model Formula",
    Predictors = "Predictors",
    logLik = "Log-Likelihood",
    AIC = "AIC",
    AICc = "AICc"
  ) %>%
  fmt_number(
    columns = c(logLik, AIC, AICc),
    decimals = 2
  ) %>%
  data_color(
    columns = AICc,
    colors = scales::col_numeric(
      palette = c("#1b7837", "#f7f7f7", "#762a83"),  # Green (best) to White to Purple (worst)
      domain = range(model_results_sorted$AICc, na.rm = TRUE)
    )
  ) %>%
  data_color(
    columns = AIC,
    colors = scales::col_numeric(
      palette = c("#1b7837", "#f7f7f7", "#762a83"),  
      domain = range(model_results_sorted$AIC, na.rm = TRUE)
    )
  ) %>%
  tab_options(
    table.font.size = px(14),
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(14),
    column_labels.font.size = px(14),
    row_group.font.size = px(14)
  )

# Display the table
model_results_table









##################################################################################################
#Now that we know phylogeny matters lets test for whether species matters in addition  
##################################################################################################

#we can't include phylogeny without species so besrt we can do is compare species to no species

m_nosp<-rma.mv(yi=betanls2_asinh, V=betanlsvar_asinh, data=all,
               mods=~1,
               random = list(~1 | study_num/substudy_num),
               method="REML",test="t")

summary(m_nosp)


m_gsp<-rma.mv(yi=betanls2_asinh, V=betanlsvar_asinh, data=all,
              mods=~1,
              random = list(~1 | study_num/substudy_num,
                            ~1 | g_sp),
              R = list(g_sp = phylo_vcv),
              method="REML",test="t")

summary(m_gsp)




AIC(m_gsp, m_nosp)


#Z test to 

# Extract species variance component and standard error
species_var <- m_gsp$sigma2[3]  # Third variance component (species-level variance)
species_se <- sqrt(species_var) # Standard deviation (square root of variance)

# Compute Z-score
z_species <- species_var / species_se

# Compute two-tailed p-value
p_species <- 2 * (1 - pnorm(abs(z_species)))  

# Print results
cat("Z-score for species variance:", z_species, "\n")
cat("P-value for species variance:", p_species, "\n")


####variance components

# Extract variance components from the model
sigma2_study <- m_gsp$sigma2[1]  # Variance among studies (expected to be ~0)
sigma2_within <- m_gsp$sigma2[2]  # Variance within species (substudy-level)
sigma2_among <- m_gsp$sigma2[3]  # Variance among species (species-level)

# Compute total variance (sum of within and among species variance)
total_variance <- sigma2_within + sigma2_among

# Compute proportion of variance at each level
prop_within <- sigma2_within / total_variance
prop_among <- sigma2_among / total_variance

# Print results
cat("Study-Level Variance: σ²_study =", round(sigma2_study, 5), "(Negligible)\n")
cat("Within-Species Variance: σ²_within =", round(sigma2_within, 3), "\n")
cat("Among-Species Variance: σ²_among =", round(sigma2_among, 3), "\n")
cat("Total Variance:", round(total_variance, 3), "\n")
cat("Proportion of Variance Within Species:", round(prop_within * 100, 1), "%\n")
cat("Proportion of Variance Among Species:", round(prop_among * 100, 1), "%\n")


# ─────────────────────────────────────────────────────────────────────────────
# Additional summary: within vs. among species variance interpretation
library(dplyr)

# Filter your data for the two focal species
wrasse_data <- all_dat2 %>%
  filter(g_sp %in% c("Thalassoma_bifasciatum", "Halichoeres_garnoti")) %>%
  select(g_sp, study_num, substudy_num, betanls2_raw_cm)

# 1. Compute species means (i.e., among-species signal)
species_means <- wrasse_data %>%
  group_by(g_sp) %>%
  summarise(species_mean = mean(betanls2_raw_cm, na.rm = TRUE))

# 2. Merge species means back into the main dataframe
wrasse_data <- wrasse_data %>%
  left_join(species_means, by = "g_sp") %>%
  mutate(deviation = betanls2_raw_cm - species_mean)

# 3. Variance components
within_species_var <- var(wrasse_data$deviation, na.rm = TRUE)
among_species_var <- var(species_means$species_mean, na.rm = TRUE)
total_var <- var(wrasse_data$betanls2_raw_cm, na.rm = TRUE)

# 4. Proportion of total variance
prop_within <- within_species_var / total_var
prop_among <- among_species_var / total_var

# 5. Output results
cat("Within-Species Variance:", round(within_species_var, 2), "\n")
cat("Among-Species Variance:", round(among_species_var, 2), "\n")
cat("Total Variance:", round(total_var, 2), "\n")
cat("Proportion Within-Species:", round(prop_within * 100, 1), "%\n")
cat("Proportion Among-Species:", round(prop_among * 100, 1), "%\n")


# ─────────────────────────────────────────────────────────────────────────────
# Example species: Thalassoma bifasciatum and Halichoeres garnoti
# ─────────────────────────────────────────────────────────────────────────────

# Filter data for these two species
wrasse_data <- all_dat2 %>%
  filter(g_sp %in% c("Thalassoma_bifasciatum", "Halichoeres_garnoti")) %>%
  select(g_sp, study_num, substudy_num, betanls2_raw_cm, betanls2_asinh) %>%
  arrange(g_sp, study_num, substudy_num)

# Summary by species
wrasse_summary <- wrasse_data %>%
  group_by(g_sp) %>%
  summarise(
    min_beta = min(betanls2_raw_cm, na.rm = TRUE),
    max_beta = max(betanls2_raw_cm, na.rm = TRUE),
    n_substudies = n()
  )

print(wrasse_summary)

# Commentary for results section
cat("\nAcross", wrasse_summary$n_substudies[wrasse_summary$g_sp == "Thalassoma_bifasciatum"],
    "substudies, *Thalassoma bifasciatum* (bluehead wrasse) showed both strongly negative and positive estimates of density dependence, ranging from",
    round(wrasse_summary$min_beta[wrasse_summary$g_sp == "Thalassoma_bifasciatum"], 2), "to",
    round(wrasse_summary$max_beta[wrasse_summary$g_sp == "Thalassoma_bifasciatum"], 2), ".\n")

cat("*Halichoeres garnoti* (yellowhead wrasse) also exhibited wide variation in beta estimates, differing by more than an order of magnitude across substudies.\n")


##################################################################################################
#Look at model output for best fit model and summarize 
##################################################################################################

all$expt_obs<-as.factor(all$expt_obs)


###Best model is just additive with no interactions 

m_all<-rma.mv(yi=betanls2_asinh, V=betanlsvar_asinh, data=all,
              mods=~logmeandensity+duration+max_len_cm+
                size_start+expt_obs,
              random = list(~1 | study_num/substudy_num, ~1 | g_sp),
              R = list(g_sp = phylo_vcv),
              method="REML",test="t")

summary(m_all)


# Extract coefficients and related statistics from the model
model_results <- data.frame(
  Term = rownames(m_all$b), # Extract predictor names
  Estimate = sinh(as.numeric(m_all$b)), # Back-transform estimates
  SE = sinh(as.numeric(m_all$se)), # Back-transform standard errors
  t_Value = as.numeric(m_all$zval), # t-values
  df = as.numeric(m_all$ddf), # Degrees of freedom
  P_Value = as.numeric(m_all$pval), # P-values
  CI = paste0("[", round(sinh(m_all$ci.lb), 3), ", ", round(sinh(m_all$ci.ub), 3), "]") # Back-transform 95% CI
)

# Convert p-values to character for formatting
model_results$P_Value <- ifelse(model_results$P_Value < 0.001, "<0.001", 
                                ifelse(model_results$P_Value < 0.01, round(model_results$P_Value, 3),
                                       ifelse(model_results$P_Value < 0.05, round(model_results$P_Value, 3),
                                              ifelse(model_results$P_Value < 0.1, round(model_results$P_Value, 3),
                                                     round(model_results$P_Value, 3)))))

# Create a polished `gt` table with back-transformed values
model_results %>%
  gt() %>%
  tab_header(
    title = md("Best-Fitting Multivariate Meta-Analysis Model"),
    subtitle = "Modeling Effect Size with Various Predictors (Back-Transformed Estimates)"
  ) %>%
  fmt_number(columns = c("Estimate", "SE", "t_Value"), decimals = 3) %>%
  cols_label(
    Term = "Predictor",
    Estimate = "Estimate (Back-Transformed)",
    SE = "Standard Error (Back-Transformed)",
    t_Value = "t-value",
    df = "Degrees of Freedom",
    P_Value = "p-value",
    CI = "95% Confidence Interval (Back-Transformed)"
  ) %>%
  tab_options(
    table.font.size = px(14),
    table.border.top.style = "solid",
    table.border.top.color = "black",
    table.border.top.width = px(2),
    table.border.bottom.style = "solid",
    table.border.bottom.color = "black",
    table.border.bottom.width = px(2),
    column_labels.font.weight = "bold",
    heading.align = "center"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = model_results$P_Value < 0.05) # Bold significant predictors
  )



##################################################################################################
#Density versus beta plot
##################################################################################################
# ─────────────────────────────────────────────────────────────────────────────
# Density vs β: Prediction Ribbon + Line Over Raw Points 
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)

# 1. Compute held-constant means
mean_duration   <- mean(all$duration,      na.rm = TRUE)
mean_size_start <- mean(all$size_start,    na.rm = TRUE)
mean_max_len_cm <- mean(all$max_len_cm,    na.rm = TRUE)

# 2. Sequence over the focal predictor
logmeandensity_vals <- seq(
  min(all$logmeandensity, na.rm = TRUE),
  max(all$logmeandensity, na.rm = TRUE),
  length.out = 100
)

# 3. Build the prediction grid (must include every term in your model formula)
pred_grid <- expand.grid(
  logmeandensity = logmeandensity_vals,
  duration       = mean_duration,
  size_start     = mean_size_start,
  max_len_cm     = mean_max_len_cm,
  expt_obs       = "Obs"               # factor level
) %>%
  # convert to the same dummy coding used in the model:
  mutate(
    expt_obs = ifelse(expt_obs == "Obs", 1, 0),
    mean_density = exp(logmeandensity)  # for plotting on the x-axis
  )

# 4. Build the newmods matrix in the exact order of your model’s moderators:
newmods <- as.matrix(pred_grid[, c(
  "logmeandensity",
  "duration",
  "max_len_cm",
  "size_start",
  "expt_obs"
)])

# 5. Predict on the arcsinh scale
pred <- predict(m_all, newmods = newmods, level = 0)

# 6. Attach back-transformed predictions & CIs
pred_grid <- pred_grid %>%
  mutate(
    pred_sinh = pred$pred,
    lb_sinh   = pred$ci.lb,
    ub_sinh   = pred$ci.ub,
    pred_beta = sinh(pred_sinh),
    lb_beta   = sinh(lb_sinh),
    ub_beta   = sinh(ub_sinh)
  )

# 7. Plot everything
steel_blue <- "#4682B4"
log_ticks  <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)

ggplot() +
  
  # raw data points
  geom_point(
    data    = all,
    mapping = aes(x = mean_density, y = betanls2_raw_cm),
    color   = steel_blue, alpha = 0.7, size = 3
  ) +
  
  # 95% CI ribbon
  geom_ribbon(
    data        = pred_grid,
    mapping     = aes(x = mean_density, ymin = lb_beta, ymax = ub_beta),
    inherit.aes = FALSE,
    fill        = "black", alpha = 0.2
  ) +
  
  # predicted β curve
  geom_line(
    data        = pred_grid,
    mapping     = aes(x = mean_density, y = pred_beta),
    inherit.aes = FALSE,
    color       = "black", size = 1.2
  ) +
  
  # zero reference
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  # scales
  scale_x_continuous(
    trans  = "log10",
    breaks = log_ticks,
    labels = scales::comma
  ) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  
  # labels & theme
  labs(
    x = "Mean Fish Density (fish/cm², log scale)",
    y = expression("Density-dependent mortality ("*beta*")")
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title      = element_text(face = "bold"),
    plot.margin     = margin(10, 10, 10, 10),
    legend.position = "none"
  )

# 8. Save the figure
ggsave(
  filename = "figures/mean_density_vs_beta_predicted_log_transformed.png",
  width    = 5,
  height   = 5,
  units    = "in",
  dpi      = 300
)









# ─────────────────────────────────────────────────────────────────────────────
# Duration vs β: Prediction Ribbon + Line Over Raw Points
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)

# 1. Held‐constant means for the other predictors
mean_logdens <- mean(all$logmeandensity, na.rm = TRUE)
mean_size    <- mean(all$size_start,      na.rm = TRUE)
mean_maxlen  <- mean(all$max_len_cm,      na.rm = TRUE)

# 2. Sequence over duration
dur_vals <- seq(
  min(all$duration, na.rm = TRUE),
  max(all$duration, na.rm = TRUE),
  length.out = 100
)

# 3. Build the prediction grid (must include every term in the model)
pred_grid_dur <- expand.grid(
  logmeandensity = mean_logdens,
  duration       = dur_vals,
  size_start     = mean_size,
  max_len_cm     = mean_maxlen,
  expt_obs       = "Obs"           # pick one level
) %>%
  # turn that into the same numeric dummy your model used:
  mutate(expt_obs = ifelse(expt_obs == "Obs", 1, 0))

# 4. Assemble newmods in the exact order of your model’s mods:
newmods_dur <- as.matrix(pred_grid_dur[, c(
  "logmeandensity",
  "duration",
  "max_len_cm",
  "size_start",
  "expt_obs"
)])

# 5. Predict on the arcsinh scale
pred_dur <- predict(m_all, newmods = newmods_dur, level = 0)

# 6. Back‐transform the fit and its 95% CI
pred_grid_dur <- pred_grid_dur %>%
  mutate(
    pred_sinh = pred_dur$pred,
    lb_sinh   = pred_dur$ci.lb,
    ub_sinh   = pred_dur$ci.ub,
    pred_beta = sinh(pred_sinh),
    lb_beta   = sinh(lb_sinh),
    ub_beta   = sinh(ub_sinh)
  )

# 7. Plot
steel_blue <- "#4682B4"

ggplot() +
  # raw data
  geom_point(
    data    = all,
    aes(x = duration, y = betanls2_raw_cm),
    color   = steel_blue, alpha = 0.7, size = 3
  ) +
  # ribbon for 95% CI
  geom_ribbon(
    data        = pred_grid_dur,
    aes(x = duration, ymin = lb_beta, ymax = ub_beta),
    inherit.aes = FALSE,
    fill        = "black", alpha = 0.2
  ) +
  # prediction line
  geom_line(
    data        = pred_grid_dur,
    aes(x = duration, y = pred_beta),
    inherit.aes = FALSE,
    color       = "black", size = 1.2
  ) +
  # zero reference
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # y-axis (asinh)
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  # labels & theme
  labs(
    x = "Study Duration (days)",
    y = expression("Density-dependent mortality ("*beta*")")
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title      = element_text(face = "bold"),
    legend.position = "none",
    plot.margin     = margin(10, 10, 10, 10)
  )

# 8. Save
ggsave(
  filename = "figures/duration_vs_beta_predicted.png",
  width    = 5,
  height   = 5,
  units    = "in",
  dpi      = 300
)




##################################################################################################
# Experimental vs. Observational: Raw + Meta‐Analytic Predictions
##################################################################################################

library(dplyr)
library(ggplot2)
library(metafor)

# 1. Compute grand means for all other covariates
mean_logdens <- mean(all$logmeandensity, na.rm = TRUE)
mean_dur     <- mean(all$duration,       na.rm = TRUE)
mean_size    <- mean(all$size_start,     na.rm = TRUE)
mean_maxlen  <- mean(all$max_len_cm,     na.rm = TRUE)

# 2. Build prediction grid with both levels of expt_obs
pred_df_eo <- expand.grid(
  expt_obs       = factor(c("Exp","Obs"), levels = c("Exp","Obs")),
  logmeandensity = mean_logdens,
  duration       = mean_dur,
  size_start     = mean_size,
  max_len_cm     = mean_maxlen
)

# 3. Construct the same moderator matrix your model saw (drop intercept)
X_new_eo <- model.matrix(
  ~ logmeandensity + duration + max_len_cm + size_start + expt_obs,
  data = pred_df_eo
)[, -1]

# 4. Predict on the arcsinh scale
pred_eo <- predict(m_all, newmods = X_new_eo, level = 0)

# 5. Back‐transform and bind into your grid
pred_df_eo <- pred_df_eo %>%
  mutate(
    fit_sinh = pred_eo$pred,
    lb_sinh  = pred_eo$ci.lb,
    ub_sinh  = pred_eo$ci.ub,
    fit_beta = sinh(fit_sinh),
    lb_beta  = sinh(lb_sinh),
    ub_beta  = sinh(ub_sinh)
  )

# 6. Colors
cols <- c(Exp = "#1E90FF", Obs = "#FF4500")

# 7. Plot
ggplot(all, aes(x = expt_obs, y = betanls2_raw_cm, fill = expt_obs)) +
  
  # raw distribution
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", alpha = 0.5) +
  geom_jitter(aes(color = expt_obs), width = 0.2, size = 3, alpha = 0.6) +
  
  # meta‐analytic 95% CI ribbon
  geom_errorbar(
    data        = pred_df_eo,
    inherit.aes = FALSE,
    aes(x = expt_obs, ymin = lb_beta, ymax = ub_beta),
    width = 0.1, size = 1.2, color = "black"
  ) +
  
  # meta‐analytic mean points
  geom_point(
    data        = pred_df_eo,
    inherit.aes = FALSE,
    aes(x = expt_obs, y = fit_beta),
    shape = 21, size = 5, fill = "black", color = "black"
  ) +
  
  # zero line & scales
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  
  # labels & theme
  labs(
    x = "Study Type",
    y = expression(beta~": strength of density-dependent mortality")
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title      = element_text(face = "bold"),
    axis.text       = element_text(color = "black"),
    plot.margin     = margin(10, 10, 10, 10)
  )

# 8. Save
ggsave(
  "figures/expt_vs_obs_beta_boxplot.png",
  width = 5, height = 5, units = "in", dpi = 300
)

# 9. Inspect your predictions
print(pred_df_eo)











library(dplyr)
library(ggplot2)
library(scales)   # for comma()

# 1) Compute means of the other moderators:
mean_logdens <- mean(all$logmeandensity, na.rm = TRUE)
mean_dur     <- mean(all$duration,       na.rm = TRUE)
mean_size    <- mean(all$size_start,     na.rm = TRUE)
mean_maxlen  <- mean(all$max_len_cm,     na.rm = TRUE)

# 2) Build Exp vs Obs “newmods” grid at those means:
pred_grid <- expand.grid(
  expt_obs       = factor(c("Exp","Obs"), levels = c("Exp","Obs")),
  logmeandensity = mean_logdens,
  duration       = mean_dur,
  size_start     = mean_size,
  max_len_cm     = mean_maxlen
)

# 3) Construct model matrix and drop intercept
X <- model.matrix(~ logmeandensity + duration + max_len_cm + size_start + expt_obs, data = pred_grid)
newmods <- X[, -1, drop = FALSE]

# 4) Predict on arcsinh scale and back-transform to β
pred_obj <- predict(m_all, newmods = newmods)

# 5) Back-transform predictions into original β units
pred_df <- pred_grid %>%
  mutate(
    pred_sinh = pred_obj$pred,
    lb_sinh   = pred_obj$ci.lb,
    ub_sinh   = pred_obj$ci.ub,
    pred_beta = sinh(pred_sinh),
    lb_beta   = sinh(lb_sinh),
    ub_beta   = sinh(ub_sinh)
  )

# 6) Join expt_obs_pairs from all_dat2
all <- all %>%
  left_join(all_dat2 %>% select(study_num, substudy_num, expt_obs_pairs),
            by = c("study_num", "substudy_num"))

# 7) Split paired vs unpaired based on expt_obs_pairs
expobs <- all %>%
  filter(expt_obs %in% c("Exp", "Obs"))

paired_expobs <- expobs %>%
  filter(!is.na(expt_obs_pairs.x))

unpaired_expobs <- expobs %>%
  filter(is.na(expt_obs_pairs.x))

# 8) Create a species-based color palette
species_cols <- setNames(
  pal2[seq_along(unique(paired_expobs$g_sp))],
  unique(paired_expobs$g_sp)
)

ggplot() +
  geom_jitter(
    data    = unpaired_expobs,
    aes(x = expt_obs, y = sinh(betanls2_asinh)),
    width   = 0.1,
    size    = 2,
    colour  = "grey70",
    alpha   = 0.5
  ) +
  geom_path(
    data   = paired_expobs,
    aes(x = expt_obs, y = sinh(betanls2_asinh),
        group = expt_obs_pairs.x,   # <- FIXED NAME
        colour = g_sp),
    size   = 0.8,
    alpha  = 0.8
  ) +
  geom_point(
    data   = paired_expobs,
    aes(x = expt_obs, y = sinh(betanls2_asinh),
        fill = g_sp,
        colour = g_sp),
    shape  = 21,
    size   = 3,
    stroke = 0.8
  ) +
  geom_errorbar(
    data   = pred_df,
    aes(x = expt_obs, ymin = lb_beta, ymax = ub_beta),
    width  = 0.2,
    size   = 1.1,
    colour = "black"
  ) +
  geom_point(
    data   = pred_df,
    aes(x = expt_obs, y = pred_beta),
    shape  = 23,
    size   = 4,
    fill   = "black",
    colour = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(name = "Species", values = species_cols) +
  scale_fill_manual(name = "Species", values = species_cols) +
  scale_y_continuous(
    name   = "Density-Dependent Mortality (β)",
    trans  = "asinh",
    breaks = c(-100, -10, -1, 0, 1, 10, 100),
    labels = comma
  ) +
  scale_x_discrete(name = "Study Type", labels = c("Exp", "Obs")) +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "right",
    axis.title       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )

##################################################################################################
#Size Start effects on beta prediction from m_all model  
##################################################################################################

# Step 1: Create a sequence of size_start values for predictions
size_start_values <- seq(
  min(all$size_start, na.rm = TRUE), 
  max(all$size_start, na.rm = TRUE), 
  length.out = 100
)

# Step 2: Create a new dataset for prediction
new_pred_data_size <- expand.grid(
  size_start      = size_start_values,
  logmeandensity  = mean(all$logmeandensity, na.rm = TRUE),
  duration        = mean(all$duration,       na.rm = TRUE),
  max_len_cm      = mean(all$max_len_cm,     na.rm = TRUE),
  expt_obs        = "Obs"  # Keep study type constant
)

# Step 3: Match factor levels for expt_obs
new_pred_data_size$expt_obs <- factor(new_pred_data_size$expt_obs, levels = levels(all$expt_obs))

# Step 4: Construct model matrix for all moderators used in m_all
model_matrix_size <- model.matrix(
  ~ logmeandensity + duration + max_len_cm + size_start + expt_obs,
  data = new_pred_data_size
)

# Step 5: Drop the intercept column
newmods_matrix_size <- model_matrix_size[, -1, drop = FALSE]

# Step 6: Predict on the arcsinh scale
predicted_values_size <- predict(m_all, newmods = newmods_matrix_size)

# Step 7: Combine predictions with the input data
pred_grid_size <- new_pred_data_size %>%
  mutate(
    pred_sinh = predicted_values_size$pred,
    lb_sinh   = predicted_values_size$ci.lb,
    ub_sinh   = predicted_values_size$ci.ub,
    pred_beta = sinh(pred_sinh),
    lb_beta   = sinh(lb_sinh),
    ub_beta   = sinh(ub_sinh)
  )

# Step 8: Plot
ggplot(pred_grid_size, aes(x = size_start, y = pred_beta)) +
  geom_ribbon(aes(ymin = lb_beta, ymax = ub_beta), fill = "black", alpha = 0.2) +
  geom_line(color = "black", size = 1.2) +
  geom_point(data = all, aes(x = size_start, y = betanls2_raw_cm),
             color = "#4682B4", alpha = 0.6, size = 3, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(
    x = "Initial Fish Size (cm)",
    y = expression("Density-dependent mortality (" * beta * ")")
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Step 9: Save the plot
ggsave(
  filename = "figures/size_start_vs_beta_predicted.png",
  width = 5,    
  height = 5,   
  units = "in",
  dpi = 300       
)







##################################################################################################
#max length and density effects with random effects 
##################################################################################################


#plot max length and log mean density regression 
# Step 1: Create a sequence of max_len_cm values for predictions
ggplot(all,aes(x=max_len_cm,y=logmeandensity))+
  geom_point()+
  theme_classic()

#fit the meta‐analytic model with both main effects & interaction:
m_int <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ max_len_cm * logmeandensity,    # includes both main effects + their interaction
  random = list(
    ~ 1 | study_num/substudy_num,
    ~ 1 | g_sp
  ),
  R      = list(g_sp = phylo_vcv),
  data   = all,
  method = "REML",
  test   = "t"
)
summary(m_int)

# 5. Extract and back‐transform the three coefficients:
coefs <- coef(summary(m_int))
results <- tibble(
  term       = rownames(coefs),
  est_asinh  = coefs[,"estimate"],
  ci_lb_asinh= coefs[,"ci.lb"],
  ci_ub_asinh= coefs[,"ci.ub"],
  
  # back‐transform via sinh()
  est         = sinh(est_asinh),
  ci_lb       = sinh(ci_lb_asinh),
  ci_ub       = sinh(ci_ub_asinh)
)

print(results)



################################################################################
# Model with max length only (no density term)
################################################################################

m_maxlen_only <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ max_len_cm,                   # only max length
  random = list(
    ~ 1 | study_num/substudy_num,
    ~ 1 | g_sp
  ),
  R      = list(g_sp = phylo_vcv),
  data   = all,
  method = "REML",
  test   = "t"
)
summary(m_maxlen_only)

# Extract and back‐transform the coefficient:
coefs_ml <- coef(summary(m_maxlen_only))
results_ml <- tibble(
  term        = rownames(coefs_ml),
  est_asinh   = coefs_ml[,"estimate"],
  ci_lb_asinh = coefs_ml[,"ci.lb"],
  ci_ub_asinh = coefs_ml[,"ci.ub"],
  
  # back‐transform via sinh()
  est          = sinh(est_asinh),
  ci_lb        = sinh(ci_lb_asinh),
  ci_ub        = sinh(ci_ub_asinh)
)

print(results_ml)


################################################################################
# Model with mean density only (no length term)
################################################################################

m_density_only <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ logmeandensity,                # only density
  random = list(
    ~ 1 | study_num/substudy_num,
    ~ 1 | g_sp
  ),
  R      = list(g_sp = phylo_vcv),
  data   = all,
  method = "REML",
  test   = "t"
)
summary(m_density_only)

# Extract and back‐transform the coefficient:
coefs_den <- coef(summary(m_density_only))
results_den <- tibble(
  term        = rownames(coefs_den),
  est_asinh   = coefs_den[,"estimate"],
  ci_lb_asinh = coefs_den[,"ci.lb"],
  ci_ub_asinh = coefs_den[,"ci.ub"],
  
  # back‐transform via sinh()
  est          = sinh(est_asinh),
  ci_lb        = sinh(ci_lb_asinh),
  ci_ub        = sinh(ci_ub_asinh)
)

print(results_den)



################################################################################
# Species effect 
################################################################################


# ──────────────────────────────────────────────────────────────────────
# 1. Prepare the data exactly as for the model (incl. study IDs!)
# ──────────────────────────────────────────────────────────────────────
all2 <- all_dat2 %>%
  filter(!is.na(betanls2_asinh),
         !is.na(betanlsvar_asinh),
         predators == "present") %>%     # correct level
  rename(beta_raw = betanls2_raw_cm) %>%
  select(g_sp, study_num, substudy_num, beta_raw, betanls2_asinh, betanlsvar_asinh) %>%
  droplevels()

# ──────────────────────────────────────────────────────────────────────
# 2. Fit the random‐effects species model (no intercept)
# ──────────────────────────────────────────────────────────────────────
m_sp_ni <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ g_sp - 1,
  random = list(~1 | study_num/substudy_num),
  data   = all2,
  method = "REML",
  test   = "t"
)

# ──────────────────────────────────────────────────────────────────────
# 3. Pull out which species actually went into that model
# ──────────────────────────────────────────────────────────────────────
model_species <- sub("^g_sp","", rownames(coef(summary(m_sp_ni))))

# ──────────────────────────────────────────────────────────────────────
# 4. Extract & back‐transform estimates + CIs
# ──────────────────────────────────────────────────────────────────────
sp_coefs <- coef(summary(m_sp_ni)) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  rename(
    est_asinh   = estimate,
    ci.lb_asinh = ci.lb,
    ci.ub_asinh = ci.ub
  ) %>%
  mutate(
    g_sp     = sub("^g_sp","",term),
    beta_est = sinh(est_asinh),
    ci.lb    = sinh(ci.lb_asinh),
    ci.ub    = sinh(ci.ub_asinh)
  ) %>%
  select(g_sp, beta_est, ci.lb, ci.ub)

# ──────────────────────────────────────────────────────────────────────
# 5. Restrict both to the modeled species
# ──────────────────────────────────────────────────────────────────────
all2     <- all2    %>% filter(g_sp %in% model_species) %>% droplevels()
sp_coefs <- sp_coefs %>% filter(g_sp %in% model_species)

# ──────────────────────────────────────────────────────────────────────
# 6. Add fish‐family for coloring
# ──────────────────────────────────────────────────────────────────────
family_df <- all_dat2 %>% distinct(g_sp, family)
all2     <- left_join(all2,     family_df, by = "g_sp")
sp_coefs <- left_join(sp_coefs, family_df, by = "g_sp")

# ──────────────────────────────────────────────────────────────────────
# 7. Order families by mean β, then species by β within each family
# ──────────────────────────────────────────────────────────────────────
family_order <- sp_coefs %>%
  group_by(family) %>%
  summarise(mean_beta = mean(beta_est)) %>%
  arrange(desc(mean_beta)) %>%
  pull(family)

species_order <- sp_coefs %>%
  mutate(family = factor(family, levels = family_order)) %>%
  arrange(family, desc(beta_est)) %>%
  pull(g_sp)

all2$g_sp     <- factor(all2$g_sp,     levels = species_order)
sp_coefs$g_sp <- factor(sp_coefs$g_sp, levels = species_order)
# ──────────────────────────────────────────────────────────────────────
# After step 7 (you already have `family_order` and `species_order`)
# ──────────────────────────────────────────────────────────────────────

# 7b) make sure `family` is a factor in the same order
all2$family    <- factor(all2$family,    levels = family_order)
sp_coefs$family <- factor(sp_coefs$family, levels = family_order)

# 8a) drop any non‐finite β’s so geom_point/geom_jitter won’t silently remove them
all2     <- all2    %>% filter(is.finite(beta_raw))
sp_coefs <- sp_coefs %>% filter(is.finite(beta_est), is.finite(ci.lb), is.finite(ci.ub))

# ──────────────────────────────────────────────────────────────────────
# After you’ve built all2, sp_coefs, family_order & species_order…
# ──────────────────────────────────────────────────────────────────────

library(ggplot2)
library(wesanderson)

# 1) Build Darjeeling1 ×9
palette_darjeeling9 <- wes_palette("Darjeeling1",
                                   type = "continuous",
                                   n    = 9)

# 2) Re-draw your plot with that palette
p_sp_darjeeling <- ggplot() +
  # raw sub‐study dots
  geom_jitter(
    data    = all2,
    aes(x = beta_raw, y = g_sp, colour = family),
    width   = 0.2, height = 0, size = 2, alpha = 0.5
  ) +
  # species‐level CIs
  geom_errorbarh(
    data    = sp_coefs,
    aes(y = g_sp, xmin = ci.lb, xmax = ci.ub, colour = family),
    height = 0, size = 0.8
  ) +
  # species‐level point estimates
  geom_point(
    data    = sp_coefs,
    aes(x = beta_est, y = g_sp, colour = family),
    shape = 17, size = 4
  ) +
  # vertical zero line
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  # asinh x‐axis
  scale_x_continuous(
    trans  = "asinh",
    breaks = c(-1000,-100,-10,-1,0,1,10,100,1000,10000),
    labels = scales::comma
  ) +
  # italic y‐axis labels (no underscores)
  scale_y_discrete(
    labels = function(x) {
      txt <- gsub("_","~", x)
      parse(text = paste0("italic(", txt, ")"))
    }
  ) +
  # Darjeeling1×9 manual palette, legend in same top→bottom order
  scale_colour_manual(
    values = palette_darjeeling9,
    breaks = family_order,
    name   = "Fish Family",
    guide  = guide_legend(reverse = TRUE)
  ) +
  labs(
    x = "Strength of Density-Dependent Mortality (β)",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y        = element_text(size = 10),
    axis.title.x       = element_text(size = 12, face = "bold"),
    legend.position    = "right",
    panel.grid.major.x = element_line(color = "gray90", linetype = "dotted")
  )

# 3) Render & save
print(p_sp_darjeeling)
ggsave("figures/species_by_family_order_beta_darjeeling9.png",
       p_sp_darjeeling,
       width = 10, height = 8,
       units = "in", dpi = 300,bg = "white")


library(dplyr)
library(stringr)
library(gt)


# Helper to clean author names (Surname, INITIALS)
fix_author_case <- function(name) {
  name <- str_trim(name)
  parts <- str_split(name, ",\\s*")[[1]]
  
  surname <- str_to_title(parts[1])
  
  initials <- if (length(parts) > 1) {
    str_replace_all(parts[2], "\\b\\w+\\b", function(x) toupper(x))
  } else {
    ""
  }
  
  str_trim(paste0(surname, if (initials != "") paste0(", ", initials)))
}


# Function to fix author name capitalization
clean_author_name <- function(name) {
  name <- trimws(name)
  parts <- unlist(str_split(name, "\\s+|,\\s*| and "))
  
  cleaned <- sapply(parts, function(word) {
    if (str_detect(word, "^[A-Z]{2,}$")) {
      word  # preserve all-uppercase initials
    } else {
      str_to_sentence(word)  # lowercase with capital first letter
    }
  })
  
  # Recombine with original comma placements
  cleaned_name <- str_replace_all(name, "\\b[A-Z]{2,}\\b", function(x) x)  # preserve initials
  cleaned_name <- str_replace_all(cleaned_name, "\\b[A-Za-z]+\\b", function(x) {
    if (str_detect(x, "^[A-Z]{2,}$")) x else str_to_sentence(x)
  })
  
  cleaned_name
}
# Construct cleaned citation and table
meta_table <- all_dat2 %>%
  mutate(
    Author = Authors,
    Year   = Publication.Year,
    Citation = ifelse(
      !is.na(Authors) & !is.na(`Article.Title`) & !is.na(`Source.Title`) & !is.na(Publication.Year),
      paste0(
        str_to_title(Authors), " (", Year, "). ",
        str_to_sentence(`Article.Title`), ". ",
        "*", str_to_title(`Source.Title`), "*."
      ),
      NA_character_
    )
  ) %>%
  select(
    Study      = study_num,
    Substudy   = substudy_num,
    Description= data.source..or.reason.for.exclusion.,
    Beta       = betanls2_raw_cm,
    Variance   = betanlsvar_raw_cm,
    Author,
    Year,
    Citation,
    DOI = DOI,
    g_sp,
    PairedStudy = paired_pred,
    ExptObsPair = expt_obs_pairs,
    Predators.= predators,
    Duration = duration,
    Exp_obs = expt_obs,
    SizeStart = size_start,
    MaxLen = max_len_cm,
    mean_density = mean_density,
    family = family,
  ) %>%
  distinct() %>%
  arrange(Year, Author) %>%
  mutate(
    Beta     = round(Beta, 3),
    Variance = round(Variance, 5)
  )

# Generate a GT table
meta_table %>%
  gt() %>%
  tab_header(
    title = md("**Summary of Studies Included in Meta-Analysis**"),
    subtitle = "Effect Size, Variance, and Metadata by Substudy"
  ) %>%
  fmt_number(columns = c("Beta", "Variance"), decimals = 3) %>%
  cols_label(
    Study       = "Study #",
    Substudy    = "Substudy #",
    Beta        = "Effect Size (β)",
    Variance    = "Variance (asinh)",
    Author      = "First Author(s)",
    Year        = "Year",
    Citation    = "Full Citation",
    PairedStudy = "Paired Predatory Treatment",
    ExptObsPair = "Exp/Obs Pair ID"
  ) %>%
  tab_options(
    table.font.size = px(13),
    column_labels.font.weight = "bold",
    heading.align = "center"
  )


write.csv(meta_table, "figures/meta_analysis_table.csv", row.names = FALSE, na = "")
