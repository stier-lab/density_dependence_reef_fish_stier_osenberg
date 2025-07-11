###############################################################################
# File: 2_estimate_beta_effect.R
# Purpose: Fit mixed-effects meta-analysis on effect sizes, back-transform estimates
# Author: Adrian Stier
# Date: 2025-07-09
###############################################################################

# ----------------------------------------------------------------------------
# 1. Setup & Preprocessing
# ----------------------------------------------------------------------------
source(here::here("code", "1_data_phylogeny_loading.R"))  # Load all_dat2
# Filter out incomplete cases for model fitting
model_data <- all_dat2 %>%
  dplyr::filter(
    !is.na(study_num),
    !is.na(substudy_num),
    !is.na(betanls2_asinh),
    !is.na(betanlsvar_asinh)
  )

# ----------------------------------------------------------------------------
# 2. Fit mixed-effects model
# ----------------------------------------------------------------------------
# Random effects: study nested within substudy
meta_model <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  random = list(~1 | study_num/substudy_num),
  data   = model_data,
  method = "REML",
  test   = "t"
)

# ----------------------------------------------------------------------------
# 3. Extract fixed-effect results and back-transform
# ----------------------------------------------------------------------------
coef_df <- coef(summary(meta_model)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "term") %>%
  dplyr::select(term, estimate, ci.lb, ci.ub)

# Inverse asinh transformation
inv_asinh <- function(x) sinh(x)

fixed_effect_results <- coef_df %>%
  dplyr::filter(term == "intrcpt") %>%
  dplyr::transmute(
    component     = "fixed_effect",
    estimate_asinh = estimate,
    ci_lb_asinh    = ci.lb,
    ci_ub_asinh    = ci.ub,
    estimate       = inv_asinh(estimate),
    ci_lower       = inv_asinh(ci.lb),
    ci_upper       = inv_asinh(ci.ub)
  )

# ----------------------------------------------------------------------------
# 4. Extract and back-transform variance components
# ----------------------------------------------------------------------------
var_comps <- summary(meta_model)$sigma2
names(var_comps) <- c("study", "substudy")

variance_results <- tibble::enframe(var_comps, name = "component", value = "var_asinh") %>%
  dplyr::mutate(
    sd_asinh = sqrt(var_asinh),
    estimate = sinh(sd_asinh)^2,
    ci_lower = NA_real_,
    ci_upper = NA_real_
  ) %>%
  dplyr::select(component, estimate, ci_lower, ci_upper)

# ----------------------------------------------------------------------------
# 5. Compile & Save Results
# ----------------------------------------------------------------------------
results_summary <- dplyr::bind_rows(fixed_effect_results, variance_results)

# Display summary
print(results_summary)

# Optional: save to CSV
# readr::write_csv(results_summary, here::here("results", "beta_meta_summary.csv"))




# ----------------------------------------------------------------------------
# 6. Combined Orchard & Beverton-Holt Figure
# ----------------------------------------------------------------------------

# Extract the model's coefficient (beta estimate)
coef_estimate <- meta_model$b

# Extract the lower and upper confidence intervals
conf_lower <- meta_model$ci.lb
conf_upper <- meta_model$ci.ub

# Back-transform the coefficient and confidence intervals from asinh transformation
back_transformed_coef <- sinh(coef_estimate) 
back_transformed_conf_lower <- sinh(conf_lower) 
back_transformed_conf_upper <- sinh(conf_upper) 



# Load necessary packages
library(ggplot2)
library(ggbeeswarm)
library(ggtext)
library(dplyr)

beta_all_all<-all_dat2%>%
  select(beta_hat,
         beta_variance_nls2,
         study_num,
         substudy_num,
         betanlsvar_raw_cm,
         betanls2_raw_cm,
         betanls2_asinh,
         betanlsvar_asinh)



#  calculate precision as 1 / variance (since you are dealing with variance)
beta_all_all <- beta_all_all %>%
  mutate(precision = 1 / betanlsvar_raw_cm)

# Define custom breaks and labels based on quantiles or a specific range
precision_breaks <- c(1e-12, 1e4, 1e6, 1e8, 1e12)
precision_labels <- c("1e-12", "1e4", "1e6", "1e8", "1e12")


######
#alternative to visualizing the plots
#####

# Load data as provided in the original script
all_data <- read.csv("data/all_studies_looped-2024-09-11.csv")
combined_results <- read.csv("output/combined_results_2024-09-17.csv")
covariates <- read.csv("output/merged-covariates-10-21-24.csv")

# Define the Beverton-Holt recruitment function
bh_f <- function(alpha, beta, t_days, x) {
  predicted_recruits <- (x * (exp(-alpha * t_days))) / (1 + (beta * x * (1 - exp(-alpha * t_days)) / alpha))
  return(predicted_recruits)
}


# Filter and merge data
usedf <- covariates %>%
  filter(use_2024 == "yes") %>%
  select(substudy_num, family,predators)%>%
  filter(predators=="present")


# Calculate quantiles for `beta` from combined_results
beta_quantiles <- quantile(combined_results$beta, probs = c(0.10, 0.5, 0.75, 0.95), na.rm = TRUE)

# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)  # Needed to render markdown/HTML in facet labels

# Filter all_data for substudies in usedf and merge with family
all_data <- inner_join(all_data, usedf, by = "substudy_num")

# Merge all_data with the combined results (alpha, beta values)
merged_data <- inner_join(all_data, combined_results, by = "substudy_num")

# Filter substudies based on the total number of replicates
filtered_data <- merged_data %>%
  group_by(substudy_num) %>%
  filter(n() >= 10) %>%  # Keep substudies with 10 or more replicates
  ungroup() %>%
  mutate(beta_strength = beta)  # Keep beta for ordering

# Step 5: Find closest beta values to the quantiles
closest_values <- sapply(beta_quantiles, function(x) {
  filtered_data$beta[which.min(abs(filtered_data$beta - x))]
})

# Step 6: Identify the substudies matching the closest beta values
matching_substudies <- filtered_data %>%
  filter(beta %in% closest_values) %>%
  distinct(substudy_num)

# Step 7: Final data with matching substudies
final_data <- filtered_data %>%
  filter(substudy_num %in% matching_substudies$substudy_num) %>%
  arrange(beta)

# Ensure final_data contains necessary variables
# Specify the desired order of substudy numbers
desired_order <- c(212, 249, 76, 256)

# Reorder `substudy_num` as a factor with the specified levels
final_data <- final_data %>%
  filter(substudy_num %in% desired_order) %>%  # Ensure only the desired substudies are included
  mutate(substudy_num = factor(substudy_num, levels = desired_order))

predicted_data <- final_data %>%
  group_by(substudy_num) %>%
  mutate(settler_range = list(seq(0, max(n0_m2, na.rm = TRUE), length.out = 100))) %>%
  unnest(cols = c(settler_range)) %>%
  mutate(predicted_recruits = bh_f(alpha, beta, t, settler_range))

# Create dynamic facet labels with **italicized species names** and beta values
facet_labels <- final_data %>%
  mutate(genus_species = gsub("_", " ", genus_species),  # Replace underscores with spaces
         facet_label = paste0("*", genus_species, "*\nBeta: ", round(beta * 1e4, 0))) %>%
  distinct(substudy_num, facet_label)

# Plot the data and predicted curves
library(ggplot2)
library(ggtext)  # Enables markdown for italic text

beverton_holt_plot <- ggplot() +
  # Scatter points with improved transparency and color
  geom_point(data = final_data, aes(x = n0_m2, y = nt_m2), 
             size = 3, alpha = 0.8, color = "#2C3E50",  # Dark blue-gray for visibility
             position = position_jitter(width = 0.05, height = 0)) +
  
  # Predicted recruitment curve with a smooth black line
  geom_line(data = predicted_data, aes(x = settler_range, y = predicted_recruits), 
            size = 1.5, color = "#E74C3C") +  # Deep red for contrast
  
  # Dashed 1:1 line in lighter gray for reference
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#95A5A6", size = 1) +
  
  # Facets with italicized species names
  facet_wrap(~substudy_num, ncol = 1, scales = "free", 
             labeller = labeller(substudy_num = setNames(facet_labels$facet_label, facet_labels$substudy_num))) +
  
  # Labels with enhanced readability
  labs(
    x = "Initial Density (Settlers)",  # You can add units if relevant, e.g., "Settlers (per m²)"
    y = "Final Density (Recruits)"
    # title = "Beverton-Holt Recruitment Model Across Species"
  ) +
  
  # Apply an elegant theme
  theme_minimal(base_size = 14) +  # Clean, modern look with larger text
  theme(
    strip.text = element_markdown(size = 14, face = "bold"),  # Italicized species names in facet labels
    axis.text = element_text(size = 12, color = "#34495E"),  # Darker gray axis text for better readability
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "#2C3E50"),  # Centered title in dark blue-gray
    panel.grid.major = element_line(color = "gray90"),  # Subtle grid for readability
    panel.grid.minor = element_line(color = "gray95"),  # Lighter minor grid lines
    panel.background = element_rect(fill = "white")  # Clean white background
  )

# Print the improved plot
print(beverton_holt_plot)





# Define custom breaks and labels based on quantiles or a specific range
precision_breaks <- c(1e-12, 1e4, 1e6, 1e8, 1e12)
precision_labels <- c("1e-12", "1e4", "1e6", "1e8", "1e12")

# Create a single-row data frame for the global estimate
global_estimate <- data.frame(
  back_transformed_coef = back_transformed_coef,
  y_position = "Beta"
)

#unique shapes for focal substudies
shape_mapping <- c("212", "249", "76", "256")  # Substudies to highlight

# Add a column for the shape
beta_all_all <- beta_all_all %>%
  mutate(
    point_shape = ifelse(as.character(substudy_num) %in% shape_mapping, 22, 16)  # 22 for matched, 16 otherwise
  )

# Define custom breaks and labels for precision scaling
precision_breaks <- c(1e-12, 1e4, 1e6, 1e8, 1e12)
precision_labels <- c(expression(10^-12), expression(10^4), expression(10^6), expression(10^8), expression(10^12))

# Define global estimate data frame
global_estimate <- data.frame(
  back_transformed_coef = back_transformed_coef,
  y_position = "Beta"
)

# Define unique shapes for focal substudies
shape_mapping <- c("212", "249", "76", "256")

# Add a column for shape mapping
beta_all_all <- beta_all_all %>%
  mutate(
    point_shape = ifelse(as.character(substudy_num) %in% shape_mapping, 22, 16)  # 22 for focal, 16 for others
  )

# Create the plot
orchard_manual <- ggplot(beta_all_all, aes(x = betanls2_raw_cm, y = "Beta", size = precision)) +
  
  # Scatter plot with jitter for better visualization
  geom_quasirandom(aes(alpha = precision), 
                   shape = 21, fill = "#377EB8", color = "#377EB8", 
                   width = 0.15, alpha = 0.7) +
  
  # Add a horizontal confidence interval line
  geom_segment(aes(x = back_transformed_conf_lower, 
                   xend = back_transformed_conf_upper, 
                   y = "Beta", yend = "Beta"), 
               color = "black", linewidth = 1.2) +
  
  # Highlight global estimate
  geom_point(data = global_estimate,
             aes(x = back_transformed_coef, y = y_position), 
             color = "black", size = 5, shape = 21, fill = "#377EB8", stroke = 1.5) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = c(-100, -10, 0, 10, 100, 1000), linetype = "dashed", color = "gray70", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
  
  # Labels
  # Labels
  labs(
    x = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", 
        cm^2,   # cm²
        ~ fish^-1, 
        ~ day^-1,
        ")"
      )
    ),
    y = "",
    size = "Precision "
  ) +
  
  
  # Apply a clean, professional theme
  theme_minimal(base_size = 16) +
  theme(
    legend.position = c(0.85, 0.85),  # Place legend in the upper-right
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  
  # Log-scale transformation for better distribution
  scale_x_continuous(
    trans = 'asinh', 
    limits = c(-500, 10000),
    breaks = c(-500, -100, -10, -1, 0, 1, 10, 100, 1000, 2000, 10000)
  ) +
  
  # Adjust point size scaling

  # 1) size scale: break on the sqrt(precision) values
  scale_size_continuous(
    name   = "Precision",
    range  = c(3, 10),                       # your point‐size range
    breaks = sqrt(precision_breaks),         # sqrt(1e-12, 1e4, …)
    labels = precision_labels                # c("1e-12","1e4",…)
  ) +
  
  # 2) guide override so points in legend look like filled circles
  guides(
    size = guide_legend(
      override.aes = list(
        shape = 21,
        fill  = "#377EB8",
        color = "#377EB8",
        alpha = 0.7,       # match your geom_quasirandom()
        stroke = 0
      )
    )
  )+
  # Flip coordinates for better visualization
  coord_flip()

# Print the plot
print(orchard_manual)




bplot2<-plot_grid(orchard_manual,beverton_holt_plot, labels = c('A', 'B'), label_size = 12,ncol=2,rel_heights = c(1, 2))

print(bplot2)



# Save the combined plot as a PDF
ggsave(
  filename = "figures/figure1_orchard_bh_plot.pdf", 
  plot = bplot2, 
  width = 10, 
  height = 11, 
  units = "in", 
  dpi = 300
)


