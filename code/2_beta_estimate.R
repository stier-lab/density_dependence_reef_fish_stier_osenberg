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