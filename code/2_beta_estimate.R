###############################################################################
# File: 2_estimate_beta_effect.R
# Purpose:
#   1) Fit mixed-effects meta-analysis on effect sizes (asinh-transformed β)
#   2) Back-transform fixed effects & variance components correctly (delta method)
#   3) Produce an Orchard + Beverton–Holt combined figure
# Author: Adrian Stier
# Date:   2025-07-09 (rev. 2025-09-25)
###############################################################################

# ----------------------------------------------------------------------------
# 1. Setup & Preprocessing
# ----------------------------------------------------------------------------
library(dplyr)
library(tibble)
library(metafor)
library(here)

# viz pkgs used later (kept for consistency with your figure styling)
library(ggplot2)
library(tidyr)
library(ggtext)
library(ggbeeswarm)
library(cowplot)

source(here::here("code", "1_data_phylogeny_loading.R"))  # loads all_dat2

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
meta_model <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  random = list(~1 | study_num/substudy_num),
  data   = model_data,
  method = "REML",
  test   = "t"
)

# ----------------------------------------------------------------------------
# 3. Fixed-effect (intercept) results & back-transform
# ----------------------------------------------------------------------------
coef_df <- coef(summary(meta_model)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "term") %>%
  dplyr::select(term, estimate, ci.lb, ci.ub)

inv_asinh <- function(x) sinh(x)

fixed_effect_results <- coef_df %>%
  dplyr::filter(term == "intrcpt") %>%
  dplyr::transmute(
    component       = "fixed_effect",
    estimate_asinh  = estimate,
    ci_lb_asinh     = ci.lb,
    ci_ub_asinh     = ci.ub,
    estimate        = inv_asinh(estimate),
    ci_lower        = inv_asinh(ci.lb),
    ci_upper        = inv_asinh(ci.ub)
  )

# convenience values for plotting later
back_transformed_coef       <- fixed_effect_results$estimate
back_transformed_conf_lower <- fixed_effect_results$ci_lower
back_transformed_conf_upper <- fixed_effect_results$ci_upper

# ----------------------------------------------------------------------------
# 4. Variance components (delta-method back-transform)
# ----------------------------------------------------------------------------
# Var(X) ≈ [cosh(mu_Y)]^2 * Var(Y), where X = sinh(Y) and mu_Y is the intercept on asinh scale
var_comps <- as.numeric(summary(meta_model)$sigma2)
names(var_comps) <- c("study", "substudy")

# Intercept (grand mean) on the asinh scale
mu_asinh  <- as.numeric(summary(meta_model)$b[1, 1])
scale_dm  <- cosh(mu_asinh)^2

variance_results <- tibble::enframe(var_comps, name = "component", value = "var_asinh") %>%
  dplyr::mutate(
    estimate = var_asinh * scale_dm,   # back-transformed variance on original scale
    ci_lower = NA_real_,               # non-trivial; left NA unless computed explicitly
    ci_upper = NA_real_
  ) %>%
  dplyr::select(component, estimate, ci_lower, ci_upper)

# ----------------------------------------------------------------------------
# 5. Compile & (optionally) Save Results
# ----------------------------------------------------------------------------
results_summary <- dplyr::bind_rows(fixed_effect_results, variance_results)
print(results_summary)
# readr::write_csv(results_summary, here::here("results", "beta_meta_summary.csv"))

# ----------------------------------------------------------------------------
# 6. Combined Orchard & Beverton–Holt Figure
# ----------------------------------------------------------------------------

# 6.1 Beverton–Holt prep (kept your formatting/colors/labels)
all_data         <- read.csv("data/all_studies_looped-2024-09-11.csv")
combined_results <- read.csv("output/combined_results_2024-09-17.csv")
covariates       <- read.csv("output/merged-covariates-10-21-24.csv")

bh_f <- function(alpha, beta, t_days, x) {
  (x * exp(-alpha * t_days)) / (1 + (beta * x * (1 - exp(-alpha * t_days)) / alpha))
}

usedf <- covariates %>%
  dplyr::filter(use_2024 == "yes", predators == "present") %>%
  dplyr::select(substudy_num, family, predators)

all_data    <- dplyr::inner_join(all_data, usedf, by = "substudy_num")
merged_data <- dplyr::inner_join(all_data, combined_results, by = "substudy_num")

filtered_data <- merged_data %>%
  dplyr::group_by(substudy_num) %>%
  dplyr::filter(n() >= 10) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_strength = beta)

beta_quantiles <- stats::quantile(combined_results$beta,
                                  probs = c(0.10, 0.5, 0.75, 0.95),
                                  na.rm = TRUE)

closest_values <- sapply(beta_quantiles, function(x) {
  filtered_data$beta[which.min(abs(filtered_data$beta - x))]
})

matching_substudies <- filtered_data %>%
  dplyr::filter(beta %in% closest_values) %>%
  dplyr::distinct(substudy_num)

final_data <- filtered_data %>%
  dplyr::filter(substudy_num %in% matching_substudies$substudy_num) %>%
  dplyr::arrange(beta)

desired_order <- c(212, 249, 76, 256)
final_data <- final_data %>%
  dplyr::filter(substudy_num %in% desired_order) %>%
  dplyr::mutate(substudy_num = factor(substudy_num, levels = desired_order))

predicted_data <- final_data %>%
  dplyr::group_by(substudy_num) %>%
  dplyr::mutate(settler_range = list(seq(0, max(n0_m2, na.rm = TRUE), length.out = 100))) %>%
  tidyr::unnest(cols = c(settler_range)) %>%
  dplyr::mutate(predicted_recruits = bh_f(alpha, beta, t, settler_range))

facet_labels <- final_data %>%
  dplyr::mutate(genus_species = gsub("_", " ", genus_species)) %>%
  dplyr::mutate(
    facet_label = paste0(
      "*", genus_species, "*\n",
      " *β* = ", round(beta * 1e4, 0), " "
    )
  ) %>%
  dplyr::distinct(substudy_num, facet_label)

beverton_holt_plot <- ggplot() +
  geom_point(
    data = final_data, aes(x = n0_m2, y = nt_m2),
    size = 3, alpha = 0.8, color = "#2C3E50",
    position = position_jitter(width = 0.05, height = 0)
  ) +
  geom_line(
    data = predicted_data, aes(x = settler_range, y = predicted_recruits),
    size = 1.5, color = "#E74C3C"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#95A5A6", linewidth = 1) +
  facet_wrap(
    ~substudy_num, ncol = 1, scales = "free",
    labeller = labeller(substudy_num = setNames(facet_labels$facet_label,
                                                facet_labels$substudy_num))
  ) +
  labs(
    x = expression(paste("Initial Density (Settlers per ", m^{-2}, ")")),
    y = expression(paste("Final Density (Recruits per ", m^{-2}, ")"))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = ggtext::element_markdown(size = 14, face = "bold"),
    axis.text  = element_text(size = 12, color = "#34495E"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    panel.background = element_rect(fill = "white")
  )

# 6.2 Orchard plot (manual) — formatting retained
precision_breaks <- c(1e-12, 1e4, 1e6, 1e8, 1e12)
precision_labels <- c(expression(10^-12), expression(10^4), expression(10^6),
                      expression(10^8), expression(10^12))

global_estimate <- data.frame(
  back_transformed_coef = back_transformed_coef,
  back_transformed_conf_lower = back_transformed_conf_lower,
  back_transformed_conf_upper = back_transformed_conf_upper,
  y_position = "Beta"
)

shape_mapping <- c("212", "249", "76", "256")

beta_all_all <- all_dat2 %>%
  dplyr::select(
    beta_hat,
    beta_variance_nls2,
    study_num,
    substudy_num,
    betanlsvar_raw_cm,
    betanls2_raw_cm,
    betanls2_asinh,
    betanlsvar_asinh
  ) %>%
  dplyr::mutate(
    precision   = 1 / beta_variance_nls2,
    point_shape = ifelse(as.character(substudy_num) %in% shape_mapping, 22, 16)
  )

orchard_manual <- ggplot(beta_all_all, aes(
  x    = betanls2_raw_cm,
  y    = "Beta",
  size = sqrt(precision)
)) +
  ggbeeswarm::geom_quasirandom(
    aes(alpha = precision),
    shape = 21, fill = "#377EB8", color = "#377EB8",
    width = 0.15, alpha = 0.7
  ) +
  geom_segment(
    data = global_estimate,
    aes(x = back_transformed_conf_lower,
        xend = back_transformed_conf_upper,
        y = y_position, yend = y_position),
    color = "black", linewidth = 1.2, inherit.aes = FALSE
  ) +
  geom_point(
    data  = global_estimate,
    aes(x = back_transformed_coef, y = y_position),
    color = "black", size = 5, shape = 21,
    fill = "#377EB8", stroke = 1.5,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", linewidth = 1) +
  geom_vline(xintercept = c(-100, -10, 0, 10, 100, 1000),
             linetype = "dashed", color = "gray70", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid",
             color = "black", linewidth = 1) +
  labs(
    x    = expression(paste("Strength of density-dependent mortality, ",
                            beta, " (", cm^2, " ", fish^{-1}, " ", day^{-1}, ")")),
    y    = NULL,
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
    ),
    alpha = "none"
  ) +
  coord_flip()

# 6.3 Combine and save figure (unchanged settings)
bplot2 <- cowplot::plot_grid(
  orchard_manual, beverton_holt_plot,
  labels = c('A', 'B'), label_size = 12, ncol = 2,
  rel_heights = c(1, 2)
)
print(bplot2)

ggplot2::ggsave(
  filename = here::here("figures", "figure1_orchard_bh_plot.pdf"),
  plot     = bplot2,
  width    = 10,
  height   = 11,
  units    = "in",
  dpi      = 300,
  device   = cairo_pdf,   # keep if Cairo is available
  bg       = "white"
)




# ----------------------------------------------------------------------------
# 7. Funnel plots: raw vs. asinh; with and without high-variance studies
# ----------------------------------------------------------------------------

# Helper: generic funnel plot on the supplied scale
make_funnel <- function(df, eff_col, var_col, mu, title, xlab,
                        point_col = "#377EB8",
                        jitter_w = 0.08, jitter_h = 0.0) {
  
  df <- df %>%
    dplyr::mutate(
      effect = .data[[eff_col]],
      se     = sqrt(.data[[var_col]])
    ) %>%
    dplyr::filter(is.finite(effect), is.finite(se), se >= 0)
  
  # y-sequence for funnel boundaries
  yseq <- seq(0, max(df$se, na.rm = TRUE), length.out = 200)
  bounds <- tibble::tibble(
    se    = yseq,
    low   = mu - 1.96 * se,
    high  = mu + 1.96 * se
  )
  
  ggplot(df, aes(x = effect, y = se)) +
    # funnel boundaries
    geom_line(data = bounds, aes(x = low,  y = se),  color = "gray40", linewidth = 0.7) +
    geom_line(data = bounds, aes(x = high, y = se),  color = "gray40", linewidth = 0.7) +
    # vertical mean line
    geom_vline(xintercept = mu, linetype = "dashed", color = "black", linewidth = 0.7) +
    # points (jittered horizontally)
    geom_point(position = position_jitter(width = jitter_w, height = jitter_h),
               shape = 21, fill = point_col, color = point_col, alpha = 0.75, size = 2.8, stroke = 0.2) +
    scale_y_reverse(expand = expansion(mult = c(0.02, 0.06))) +  # standard funnel: small SE at top
    labs(x = xlab, y = "Standard error", title = title) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title        = element_text(face = "bold", size = 14),
      axis.title        = element_text(face = "bold", size = 14),
      axis.text         = element_text(size = 12, color = "#34495E"),
      panel.grid.major  = element_line(color = "gray90"),
      panel.grid.minor  = element_line(color = "gray95"),
      panel.background  = element_rect(fill = "white", color = NA)
    )
}

# 7.1 Set up data, means, and high-variance filters
mu_asinh <- as.numeric(summary(meta_model)$b[1, 1])
mu_raw   <- back_transformed_coef[1]  # from Section 3 (sinh of intercept CI mid)

# Define "high variance" as the top 5% of sampling variances on each scale
q_raw   <- stats::quantile(all_dat2$betanlsvar_raw_cm,   probs = 0.95, na.rm = TRUE)
q_asinh <- stats::quantile(all_dat2$betanlsvar_asinh,    probs = 0.95, na.rm = TRUE)

df_raw_all   <- all_dat2 %>% dplyr::filter(is.finite(betanls2_raw_cm),  is.finite(betanlsvar_raw_cm))
df_raw_trim  <- df_raw_all  %>% dplyr::filter(betanlsvar_raw_cm   <= q_raw)

df_tr_all    <- all_dat2 %>% dplyr::filter(is.finite(betanls2_asinh), is.finite(betanlsvar_asinh))
df_tr_trim   <- df_tr_all   %>% dplyr::filter(betanlsvar_asinh    <= q_asinh)

# 7.2 Build the four funnel plots
p_raw_all <- make_funnel(
  df      = df_raw_all,
  eff_col = "betanls2_raw_cm",
  var_col = "betanlsvar_raw_cm",
  mu      = mu_raw,
  title   = "Funnel (raw scale): all studies",
  xlab    = expression(paste(beta, " (", cm^2, " ", fish^{-1}, " ", day^{-1}, ")"))
)

p_raw_trim <- make_funnel(
  df      = df_raw_trim,
  eff_col = "betanls2_raw_cm",
  var_col = "betanlsvar_raw_cm",
  mu      = mu_raw,
  title   = "Funnel (raw scale): excluding top 5% variances",
  xlab    = expression(paste(beta, " (", cm^2, " ", fish^{-1}, " ", day^{-1}, ")"))
)

p_tr_all <- make_funnel(
  df      = df_tr_all,
  eff_col = "betanls2_asinh",
  var_col = "betanlsvar_asinh",
  mu      = mu_asinh,
  title   = "Funnel (asinh scale): all studies",
  xlab    = expression(paste("asinh(", beta, ")"))
)

p_tr_trim <- make_funnel(
  df      = df_tr_trim,
  eff_col = "betanls2_asinh",
  var_col = "betanlsvar_asinh",
  mu      = mu_asinh,
  title   = "Funnel (asinh scale): excluding top 5% variances",
  xlab    = expression(paste("asinh(", beta, ")"))
)

# 7.3 Arrange and save
funnel_grid <- cowplot::plot_grid(
  p_raw_all, p_raw_trim, p_tr_all, p_tr_trim,
  labels = c("A", "B", "C", "D"),
  label_size = 12, ncol = 2
)
print(funnel_grid)

ggplot2::ggsave(
  filename = here::here("figures", "figure2_funnel_plots.pdf"),
  plot     = funnel_grid,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300,
  device   = cairo_pdf,
  bg       = "white"
)

# (Optional PNG for quick looks)
# ggplot2::ggsave(
#   filename = here::here("figures", "figure2_funnel_plots.png"),
#   plot     = funnel_grid,
#   width    = 10, height = 8, units = "in", dpi = 300, bg = "white"
# )


