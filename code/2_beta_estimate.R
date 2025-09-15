# ==============================================================================
# File: 8_sensitivity_core_beta.R
# Title: Core β Diagnostics — matches 2_estimate_beta_effect.R variance usage
# Author: Stier Lab (Adrian C. Stier), with assistant support
# Date: 2025-09-15
#
# PURPOSE
#   Run reviewer-focused diagnostics for the core mixed-effects meta-analysis:
#     • Heterogeneity (tau^2, sigma^2, QE/QEp, pseudo-I^2)
#     • Leave-one-out (CSV + plot)
#     • Influence diagnostics (Cook’s D, hat, max|DFBETAS| + plot)
#     • Small-study effects: funnel plot + Egger test (trim-and-fill to SI)
#
# MODEL (identical inputs to 2_estimate_beta_effect.R)
#   yi = betanls2_asinh
#   V  = betanlsvar_asinh          # (sampling variance in the asinh space)
#   random = ~ 1 | study_num/substudy_num
#
# OUTPUTS
#   figs/sensitivity_core_beta/baseline_beta/
#       funnel.(png|pdf), influence.(png|pdf), leave_one_out.(png|pdf)
#   tables/sensitivity_core_beta/baseline_beta/
#       heterogeneity.csv, egger_test.csv, trimfill_SI.csv,
#       influence_metrics.csv, leave_one_out.csv, variance_diagnostics.csv,
#       intercept_backtransform.csv, sessionInfo_baseline_beta.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(forcats); library(broom); library(janitor)
  library(metafor); library(clubSandwich)
  library(here)
})

# ---------------------------- Paths & Helpers ---------------------------------
ROOT_FIG <- "figs/sensitivity_core_beta"
ROOT_TAB <- "tables/sensitivity_core_beta"
dir.create(ROOT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(ROOT_TAB, recursive = TRUE, showWarnings = FALSE)

mk_dirs <- function(label) {
  fig_dir <- file.path(ROOT_FIG, label); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  tab_dir <- file.path(ROOT_TAB, label); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  list(fig = fig_dir, tab = tab_dir)
}

# Publication-quality saving: PNG + vector PDF
save_plot_dual <- function(plot, outfile_base, width=7, height=5, dpi=400, bg="white") {
  ggplot2::ggsave(paste0(outfile_base, ".png"), plot, width=width, height=height, dpi=dpi, bg=bg)
  ggplot2::ggsave(paste0(outfile_base, ".pdf"), plot, width=width, height=height,
                  dpi=600, bg=bg, device = cairo_pdf)
}

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
}

# Conservative variance floor, used only for fitting/plots (raw V untouched)
compute_vi_floor <- function(v) {
  pos <- v[is.finite(v) & v > 0]
  if (!length(pos)) return(NA_real_)
  max(1e-12,
      min(min(pos, na.rm=TRUE) * 1e-3,
          stats::quantile(pos, probs = 0.005, na.rm = TRUE, names = FALSE)))
}

with_v_fit <- function(dat, vi_col = "betanlsvar_asinh") {
  v_raw <- dat[[vi_col]]
  floor_val <- compute_vi_floor(v_raw)
  dat %>% mutate(
    v_fit   = if (!is.na(floor_val)) pmax(.data[[vi_col]], floor_val) else .data[[vi_col]],
    v_floor = floor_val
  )
}

# ---------------------------- Load analysis data ------------------------------
# Must match 2_estimate_beta_effect.R
source(here::here("code", "1_data_phylogeny_loading.R"))  # loads all_dat2

model_data <- all_dat2 %>%
  filter(
    !is.na(study_num),
    !is.na(substudy_num),
    !is.na(betanls2_asinh),
    !is.na(betanlsvar_asinh)
  ) %>%
  # add v_fit used for *fitting and bias plots only*
  with_v_fit("betanlsvar_asinh")

# ----------------------------- Fit core model ---------------------------------
fit_core <- tryCatch(
  metafor::rma.mv(
    yi     = betanls2_asinh,
    V      = v_fit,                       # stays in asinh space (variance)
    random = list(~ 1 | study_num/substudy_num),
    data   = model_data,
    method = "REML",
    test   = "t"
  ),
  error = function(e) { message("Core fit failed: ", e$message); NULL }
)

# Back-transform fixed intercept (context only, not used in diagnostics)
inv_asinh <- function(x) sinh(x)
if (!is.null(fit_core)) {
  cf <- coef(summary(fit_core))
  intr <- as_tibble(cf, rownames = "term") %>%
    filter(term == "intrcpt") %>%
    transmute(
      component        = "fixed_effect",
      estimate_asinh   = estimate,
      ci_lb_asinh      = ci.lb,
      ci_ub_asinh      = ci.ub,
      estimate         = inv_asinh(estimate),
      ci_lower         = inv_asinh(ci.lb),
      ci_upper         = inv_asinh(ci.ub)
    )
  # Save for SI / narrative
  dirs <- mk_dirs("baseline_beta")
  readr::write_csv(intr, file.path(dirs$tab, "intercept_backtransform.csv"))
}

# --------------------------- Heterogeneity table ------------------------------
dirs <- mk_dirs("baseline_beta")
fig_dir <- dirs$fig; tab_dir <- dirs$tab

if (!is.null(fit_core)) {
  het <- tibble(
    k            = nrow(model_data),
    tau2_resid   = unname(fit_core$tau2),
    QE           = unname(fit_core$QE),
    QEp          = unname(fit_core$QEp),
    sigma2       = paste(round(fit_core$sigma2, 6), collapse = ";"),
    sigma2_names = paste(fit_core$s.names, collapse = ";"),
    I2_overall   = suppressWarnings(tryCatch(metafor::i2(fit_core)$I2, error = function(e) NA_real_))
  )
  readr::write_csv(het, file.path(tab_dir, "heterogeneity.csv"))
} else {
  readr::write_csv(tibble(note="Intercept fit failed; heterogeneity unavailable."),
                   file.path(tab_dir, "heterogeneity.csv"))
}

# ----------------------- Variance diagnostics (SI transparency) ----------------
vard <- tibble(
  k        = nrow(model_data),
  v_floor  = unique(model_data$v_floor),
  V_min    = min(model_data$betanlsvar_asinh, na.rm=TRUE),
  V_q05    = quantile(model_data$betanlsvar_asinh, 0.05, na.rm=TRUE),
  V_med    = median(model_data$betanlsvar_asinh, na.rm=TRUE),
  V_q95    = quantile(model_data$betanlsvar_asinh, 0.95, na.rm=TRUE),
  V_max    = max(model_data$betanlsvar_asinh, na.rm=TRUE)
)
readr::write_csv(vard, file.path(tab_dir, "variance_diagnostics.csv"))

# -------------------------- Influence diagnostics -----------------------------
if (!is.null(fit_core)) {
  cooks_cutoff <- 4 / max(5, nrow(model_data))
  
  infl <- tryCatch(metafor::influence.rma.mv(fit_core), error=function(e) NULL)
  if (!is.null(infl)) {
    inf_df <- tibble(
      idx     = seq_along(infl$cook.d),
      cooks_d = infl$cook.d,
      hat     = infl$hat,
      dfbetas = apply(infl$dfbs, 1, function(x) max(abs(x), na.rm = TRUE))
    )
    readr::write_csv(inf_df, file.path(tab_dir, "influence_metrics.csv"))
    
    p_inf <- inf_df %>%
      ggplot(aes(cooks_d, dfbetas)) +
      geom_point(size = 1.8, alpha = .9) +
      geom_vline(xintercept = cooks_cutoff, linetype = "dashed") +
      labs(x = "Cook's distance", y = "max|DFBETAS|", title = "Influence diagnostics") +
      theme_pub()
    save_plot_dual(p_inf, file.path(fig_dir, "influence"))
  }
  
  # ------------------------------ Leave-one-out -------------------------------
  loo <- tryCatch(metafor::leave1out.rma.mv(fit_core), error=function(e) NULL)
  if (!is.null(loo) && length(loo$estimate) == nrow(model_data)) {
    loo_df <- tibble(idx = seq_along(loo$estimate),
                     loo_est = as.numeric(loo$estimate))
    readr::write_csv(loo_df, file.path(tab_dir, "leave_one_out.csv"))
    
    p_loo <- loo_df %>%
      ggplot(aes(idx, loo_est)) +
      geom_point(size = 1.6, alpha = .9) +
      geom_hline(yintercept = as.numeric(fit_core$b), linetype = "dashed") +
      labs(x = "Left-out index", y = expression(beta[asinh]),
           title = "Leave-one-out β (asinh space)") +
      theme_pub()
    save_plot_dual(p_loo, file.path(fig_dir, "leave_one_out"))
  }
}

# ------------------------ Small-study effects (bias) --------------------------
# NB: use rma.uni on the SAME effect/variance used above (asinh space with v_fit)
uni <- tryCatch(
  metafor::rma.uni(yi = betanls2_asinh, vi = v_fit, method = "REML", data = model_data),
  error = function(e) NULL
)

if (!is.null(uni)) {
  # symmetric x-limits around 0 using 99th percentile of |yi|
  xlim_max <- stats::quantile(abs(model_data$betanls2_asinh), probs = 0.99, na.rm = TRUE)
  xlim_max <- max(1e-6, as.numeric(xlim_max))
  xlim_use <- c(-xlim_max, xlim_max)
  
  # base metafor funnel (it draws CI wedges correctly), save PNG + PDF
  png(file.path(fig_dir, "funnel.png"), width = 1800, height = 1400, res = 400, bg = "white")
  par(mar = c(4.2, 4.2, 3.0, 1.0))
  metafor::funnel(uni, main = "Funnel plot (asinh space)",
                  xlim = xlim_use, refline = 0, shade = c("gray95", "gray90"))
  dev.off()
  
  grDevices::cairo_pdf(file.path(fig_dir, "funnel.pdf"), width = 7, height = 5, family = "Helvetica")
  par(mar = c(4.2, 4.2, 3.0, 1.0))
  metafor::funnel(uni, main = "Funnel plot (asinh space)",
                  xlim = xlim_use, refline = 0, shade = c("gray95", "gray90"))
  dev.off()
  
  # Egger regression (small-study effects)
  egger <- tryCatch(metafor::regtest(uni, model = "lm"), error=function(e) NULL)
  if (!is.null(egger)) {
    readr::write_csv(tibble(z = unname(egger$zval), p = unname(egger$pval)),
                     file.path(tab_dir, "egger_test.csv"))
  }
  
  # Trim-and-fill (to SI)
  tf <- tryCatch(metafor::trimfill(uni), error = function(e) NULL)
  if (!is.null(tf)) readr::write_csv(broom::tidy(tf), file.path(tab_dir, "trimfill_SI.csv"))
}

# ------------------------------- Session info ---------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(tab_dir, "sessionInfo_baseline_beta.txt"))

message("Core β diagnostics complete. See: ", ROOT_FIG, " / ", ROOT_TAB)