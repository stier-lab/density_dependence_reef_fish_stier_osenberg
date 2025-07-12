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

# End of 2_estimate_beta_effect.R 


###############################################################################
# File: 2_estimate_beta_effect.R
# Purpose:
#   1. Fit mixed-effects meta-analysis on effect sizes (asinh-transformed β)
#   2. Back-transform estimates
#   3. Produce an Orchard + Beverton-Holt combined figure (Step 6)
# Author: Adrian Stier
# Date:   2025-07-09 (rev. 2025-07-10)
###############################################################################

# ----------------------------------------------------------------------------
# 1. Setup & Preprocessing
# ----------------------------------------------------------------------------
library(dplyr)
library(tibble)
library(metafor)
library(here)

source(here("code", "1_data_phylogeny_loading.R"))  # loads all_dat2

model_data <- all_dat2 %>%
  filter(
    !is.na(study_num),
    !is.na(substudy_num),
    !is.na(betanls2_asinh),
    !is.na(betanlsvar_asinh)
  )

# ----------------------------------------------------------------------------
# 2. Fit mixed-effects model
# ----------------------------------------------------------------------------
meta_model <- rma.mv(
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
  rownames_to_column("term") %>%
  select(term, estimate, ci.lb, ci.ub)

inv_asinh <- function(x) sinh(x)

fixed_effect_results <- coef_df %>%
  filter(term == "intrcpt") %>%
  transmute(
    component       = "fixed_effect",
    estimate_asinh  = estimate,
    ci_lb_asinh     = ci.lb,
    ci_ub_asinh     = ci.ub,
    estimate        = inv_asinh(estimate),
    ci_lower        = inv_asinh(ci.lb),
    ci_upper        = inv_asinh(ci.ub)
  )

# ----------------------------------------------------------------------------
# 4. Extract and back-transform variance components
# ----------------------------------------------------------------------------
var_comps <- summary(meta_model)$sigma2
names(var_comps) <- c("study", "substudy")

variance_results <- enframe(var_comps, name = "component", value = "var_asinh") %>%
  mutate(
    sd_asinh = sqrt(var_asinh),
    estimate = sinh(sd_asinh)^2,
    ci_lower = NA_real_,
    ci_upper = NA_real_
  ) %>%
  select(component, estimate, ci_lower, ci_upper)

# ----------------------------------------------------------------------------
# 5. Compile & Save Results
# ----------------------------------------------------------------------------
results_summary <- bind_rows(fixed_effect_results, variance_results)
print(results_summary)
# readr::write_csv(results_summary, here("results", "beta_meta_summary.csv"))

# ----------------------------------------------------------------------------
# 6. Combined Orchard & Beverton-Holt Figure
# ----------------------------------------------------------------------------

# Load required packages for visualization
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)
library(ggbeeswarm)
library(cowplot)

# 6.1 Load and prepare Beverton-Holt data
all_data         <- read.csv("data/all_studies_looped-2024-09-11.csv")
combined_results <- read.csv("output/combined_results_2024-09-17.csv")
covariates       <- read.csv("output/merged-covariates-10-21-24.csv")

bh_f <- function(alpha, beta, t_days, x) {
  (x * exp(-alpha * t_days)) / (1 + (beta * x * (1 - exp(-alpha * t_days)) / alpha))
}

usedf <- covariates %>%
  filter(use_2024 == "yes", predators == "present") %>%
  select(substudy_num, family, predators)

all_data    <- inner_join(all_data, usedf, by = "substudy_num")
merged_data <- inner_join(all_data, combined_results, by = "substudy_num")

filtered_data <- merged_data %>%
  group_by(substudy_num) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  mutate(beta_strength = beta)

beta_quantiles <- quantile(combined_results$beta,
                           probs = c(0.10, 0.5, 0.75, 0.95),
                           na.rm = TRUE)

closest_values <- sapply(beta_quantiles, function(x) {
  filtered_data$beta[which.min(abs(filtered_data$beta - x))]
})

matching_substudies <- filtered_data %>%
  filter(beta %in% closest_values) %>%
  distinct(substudy_num)

final_data <- filtered_data %>%
  filter(substudy_num %in% matching_substudies$substudy_num) %>%
  arrange(beta)

desired_order <- c(212, 249, 76, 256)
final_data <- final_data %>%
  filter(substudy_num %in% desired_order) %>%
  mutate(substudy_num = factor(substudy_num, levels = desired_order))

predicted_data <- final_data %>%
  group_by(substudy_num) %>%
  mutate(settler_range = list(seq(0, max(n0_m2, na.rm = TRUE), length.out = 100))) %>%
  unnest(cols = c(settler_range)) %>%
  mutate(predicted_recruits = bh_f(alpha, beta, t, settler_range))

facet_labels <- final_data %>%
  mutate(genus_species = gsub("_", " ", genus_species),
         facet_label = paste0("*", genus_species, "*\nBeta: ", round(beta * 1e4, 0))) %>%
  distinct(substudy_num, facet_label)

beverton_holt_plot <- ggplot() +
  geom_point(data = final_data, aes(x = n0_m2, y = nt_m2),
             size = 3, alpha = 0.8, color = "#2C3E50",
             position = position_jitter(width = 0.05, height = 0)) +
  geom_line(data = predicted_data, aes(x = settler_range, y = predicted_recruits),
            size = 1.5, color = "#E74C3C") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#95A5A6", size = 1) +
  facet_wrap(~substudy_num, ncol = 1, scales = "free",
             labeller = labeller(substudy_num = setNames(facet_labels$facet_label, facet_labels$substudy_num))) +
  labs(
    x = expression(paste("Initial Density (Settlers per ", m^{-2}, ")")),
    y = expression(paste("Final Density (Recruits per ", m^{-2}, ")"))
  )+
  theme_minimal(base_size = 14) +
  theme(strip.text = element_markdown(size = 14, face = "bold"),
        axis.text = element_text(size = 12, color = "#34495E"),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white"))

print(beverton_holt_plot)

# 6.2 Prepare and draw Orchard plot manually
precision_breaks <- c(1e-12, 1e4, 1e6, 1e8, 1e12)
precision_labels <- c(expression(10^-12), expression(10^4), expression(10^6), expression(10^8), expression(10^12))

global_estimate <- data.frame(
  back_transformed_coef = back_transformed_coef,
  y_position = "Beta"
)

shape_mapping <- c("212", "249", "76", "256")

beta_all_all <- beta_all_all %>%
  mutate(point_shape = ifelse(as.character(substudy_num) %in% shape_mapping, 22, 16))

orchard_manual <- ggplot(beta_all_all, aes(
  x = betanls2_raw_cm,
  y = "Beta",
  size = sqrt(precision)
)) +
  geom_quasirandom(aes(alpha = precision),
                   shape = 21, fill = "#377EB8", color = "#377EB8",
                   width = 0.15, alpha = 0.7) +
  geom_segment(aes(x = back_transformed_conf_lower,
                   xend = back_transformed_conf_upper,
                   y = "Beta", yend = "Beta"),
               color = "black", linewidth = 1.2) +
  geom_point(data = global_estimate,
             aes(x = back_transformed_coef, y = y_position),
             color = "black", size = 5, shape = 21,
             fill = "#377EB8", stroke = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", linewidth = 1) +
  geom_vline(xintercept = c(-100, -10, 0, 10, 100, 1000),
             linetype = "dashed", color = "gray70", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid",
             color = "black", linewidth = 1) +
  labs(
    x = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", cm^2, ~ fish^-1, ~ day^-1, ")"
      )
    ),
    y = NULL,
    size = "Precision (1/SE)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.y      = element_blank(),
    legend.position   = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key        = element_blank(),
    axis.text         = element_text(size = 14),
    axis.title.x      = element_text(size = 16, face = "bold"),
    panel.grid.major  = element_line(color = "gray90"),
    panel.grid.minor  = element_blank(),
    panel.background  = element_rect(fill = "white")
  ) +
  scale_x_continuous(
    trans  = "asinh",
    limits = c(-500, 10000),
    breaks = c(-500, -100, -10, -1, 0, 1, 10, 100, 1000, 2000, 10000)
  ) +
  scale_size_continuous(
    name   = "Precision",
    range  = c(3, 10),
    breaks = sqrt(precision_breaks),
    labels = precision_labels
  ) +
  guides(
    size = guide_legend(
      override.aes = list(
        shape  = 16,
        fill   = "#377EB8",
        color  = "#377EB8",
        stroke = 0
      )
    )
  ) +
  coord_flip()

print(orchard_manual)


# 6.3 Combine and save
bplot2 <- plot_grid(orchard_manual, beverton_holt_plot,
                    labels = c('A', 'B'), label_size = 12, ncol = 2,
                    rel_heights = c(1, 2))
print(bplot2)

ggsave(
  filename = here("figures", "figure1_orchard_bh_plot.pdf"),
  plot     = bplot2,
  width    = 10,
  height   = 11,
  units    = "in",
  dpi      = 300
)