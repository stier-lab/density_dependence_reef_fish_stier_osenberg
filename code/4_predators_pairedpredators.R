###############################################################################
# File: 4_predators_pairedpredators.R
# Purpose: Quantify how predator presence alters density‐dependence (β),
#          accounting for phylogeny and paired‐study design
# Author: Adrian Stier (rev. 2025‐07‐09)
###############################################################################

# ----------------------------------------------------------------------------
# 1. Dependencies & Data Load
# ----------------------------------------------------------------------------
source(here::here("code", "0_libraries.R"))               # core packages
source(here::here("code", "1_data_phylogeny_loading.R"))  # loads all_dat2, pruned_tree, etc.
source(here::here("code", "2_beta_estimate.R"))           # computes betanls2_asinh, etc.

# Create output directories if needed
dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("results"), showWarnings = FALSE, recursive = TRUE)


# ──────────────────────────────────────────────────────────────────────────────
# 1. Predator‐presence summary
# ──────────────────────────────────────────────────────────────────────────────
predator_counts <- all_dat2 %>%
  count(predators) %>%
  mutate(proportion = n / sum(n))

print(predator_counts)

cat("Fraction with predators present: ",
    round(predator_counts$proportion[predator_counts$predators == "present"], 3),
    "\n")
cat("Fraction with predators absent:  ",
    round(predator_counts$proportion[predator_counts$predators == "absent"], 3),
    "\n")


# ──────────────────────────────────────────────────────────────────────────────
# 2. Build & trim phylogenetic VCV to match data
# ──────────────────────────────────────────────────────────────────────────────
phylo_vcv_full <- ape::vcv(pruned_tree, corr = TRUE)

good_sp    <- intersect(all_dat2$g_sp, rownames(phylo_vcv_full))
model_data <- all_dat2 %>%
  filter(
    g_sp %in% good_sp,
    !is.na(betanls2_asinh),
    !is.na(betanlsvar_asinh)
  )
phylo_vcv <- phylo_vcv_full[good_sp, good_sp]

missing_sp <- setdiff(unique(all_dat2$g_sp), good_sp)
if (length(missing_sp) > 0) {
  warning("Dropped species not in tree: ",
          paste(missing_sp, collapse = ", "))
}


# ──────────────────────────────────────────────────────────────────────────────
# 3. Mixed‐effects meta‐analysis (null vs predators)
# ──────────────────────────────────────────────────────────────────────────────
m_null <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ 1,
  random = list(~1 | study_num/substudy_num,
                ~1 | g_sp),
  R      = list(g_sp = phylo_vcv),
  data   = model_data,
  method = "REML"
)

m_pred <- update(m_null, mods = ~ predators, test = "t")
print(AIC(m_null, m_pred))


# ──────────────────────────────────────────────────────────────────────────────
# 4. Extract & back‐transform predator‐model coefficients
# ──────────────────────────────────────────────────────────────────────────────
cfs        <- coef(summary(m_pred))
β_abs_log  <- cfs["intrcpt",         "estimate"]
β_pre_log  <- cfs["predatorspresent","estimate"]
ci_abs_log <- c(m_pred$ci.lb[1], m_pred$ci.ub[1])
ci_pre_log <- c(m_pred$ci.lb[2], m_pred$ci.ub[2])

β_abs  <- sinh(β_abs_log)
β_pre  <- sinh(β_abs_log + β_pre_log)
CI_abs <- sinh(ci_abs_log)
CI_pre <- sinh(β_abs_log + ci_pre_log)
pct_ch <- 100 * (β_pre - β_abs) / abs(β_abs)

cat(sprintf(
  "Meta β no predators:    %.3f [%.3f, %.3f]\n",
  β_abs, CI_abs[1], CI_abs[2]
))
cat(sprintf(
  "Meta β with predators:  %.3f [%.3f, %.3f]\n",
  β_pre, CI_pre[1], CI_pre[2]
))
cat(sprintf("Percent change:         %.1f%%\n", pct_ch))


# ──────────────────────────────────────────────────────────────────────────────
# 5. Summary table (gt)
# ──────────────────────────────────────────────────────────────────────────────
res_tbl <- tibble::tribble(
  ~Metric,                  ~Estimate,                                                       ~`95% CI`,
  "Studies w/ predators",   sprintf("%d", predator_counts$n[predator_counts$predators == "present"]),    "",
  "Studies w/o predators",  sprintf("%d", predator_counts$n[predator_counts$predators == "absent"]),     "",
  "Fraction present",       sprintf("%.3f", predator_counts$proportion[predator_counts$predators == "present"]), "",
  "Fraction absent",        sprintf("%.3f", predator_counts$proportion[predator_counts$predators == "absent"]),   "",
  "β (no predators)",       sprintf("%.3f", β_abs),                                            sprintf("[%.3f, %.3f]", CI_abs[1], CI_abs[2]),
  "β (with predators)",     sprintf("%.3f", β_pre),                                            sprintf("[%.3f, %.3f]", CI_pre[1], CI_pre[2]),
  "Pct change",             sprintf("%.1f%%", pct_ch),                                         ""
)

res_tbl %>%
  gt::gt() %>%
  gt::tab_header(
    title    = md("**Predator-Presence Meta-Analysis**"),
    subtitle = "Effect sizes & study counts"
  ) %>%
  gt::cols_label(
    Metric   = "Metric",
    Estimate = "Estimate",
    `95% CI` = "95% Confidence Interval"
  ) %>%
  gt::fmt_missing(columns = everything(), missing_text = "") %>%
  gt::tab_options(
    table.font.size           = px(14),
    heading.align             = "center",
    column_labels.font.weight = "bold"
  ) %>%
  gt::gtsave(here::here("results", "predator_presence_summary.png"))


# ──────────────────────────────────────────────────────────────────────────────
# 6. All‐vs‐paired studies comparison
# ──────────────────────────────────────────────────────────────────────────────
summarize_pred_effect <- function(mod) {
  cf    <- coef(summary(mod))
  i_log <- cf["intrcpt",         "estimate"]
  m_log <- cf["predatorspresent","estimate"]
  b0    <- sinh(i_log)
  b1    <- sinh(i_log + m_log)
  ci0   <- sinh(c(mod$ci.lb[1], mod$ci.ub[1]))
  ci1   <- sinh(i_log + c(mod$ci.lb[2], mod$ci.ub[2]))
  
  tibble(
    set      = NA_character_,
    beta0    = b0,
    ci0_lo   = ci0[1],
    ci0_hi   = ci0[2],
    beta1    = b1,
    ci1_lo   = ci1[1],
    ci1_hi   = ci1[2],
    pct_mean = 100 * (b1 - b0) / abs(b0)
  )
}

m_all      <- m_pred
paired_df  <- model_data %>% filter(paired_pred == "paired")
m_paired   <- update(m_all, data = paired_df)

sum_all    <- summarize_pred_effect(m_all)    %>% mutate(set = "All studies")
sum_paired <- summarize_pred_effect(m_paired) %>% mutate(set = "Paired only")

comp_tbl <- bind_rows(sum_all, sum_paired) %>%
  transmute(
    `Study set`   = set,
    `β no pred`   = sprintf("%.2f [%.2f, %.2f]", beta0, ci0_lo, ci0_hi),
    `β with pred` = sprintf("%.2f [%.2f, %.2f]", beta1, ci1_lo, ci1_hi),
    `% inc (mean)` = sprintf("%.1f%%", pct_mean)
  )

comp_tbl %>%
  gt::gt() %>%
  gt::tab_header(
    title    = md("**Predator‐Presence: All vs Paired**"),
    subtitle = "Meta‐analytic β estimates"
  ) %>%
  gt::cols_label(
    `Study set`   = "Study set",
    `β no pred`   = md("**β** (no predators)"),
    `β with pred` = md("**β** (with predators)")
  ) %>%
  gt::tab_options(heading.align = "center", table.font.size = px(14)) %>%
  gt::gtsave(here::here("results", "predator_all_vs_paired.png"))


# ──────────────────────────────────────────────────────────────────────────────
# 7. Boxplot + meta‐predictions
# ──────────────────────────────────────────────────────────────────────────────
new_pred <- tibble(
  predators = factor(c("absent", "present"), levels = c("absent", "present"))
)
X_new <- model.matrix(~ predators, data = new_pred)
preds <- predict(m_pred, newmods = X_new[, "predatorspresent", drop = FALSE])

new_pred <- new_pred %>%
  mutate(
    beta_hat = sinh(preds$pred),
    ci_lo    = sinh(preds$ci.lb),
    ci_hi    = sinh(preds$ci.ub)
  )

pred_cols <- c(absent = "#1E90FF", present = "#FF4500")

p1 <- ggplot(model_data, aes(predators, betanls2_raw_cm)) +
  geom_boxplot(aes(fill = predators), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = predators), width = 0.2, size = 2) +
  geom_point(
    data = new_pred,
    aes(predators, beta_hat),
    shape = 21, size = 5, color = "black", stroke = 1
  ) +
  geom_errorbar(
    data = new_pred,
    aes(predators, beta_hat, ymin = ci_lo, ymax = ci_hi),
    width = 0.1, size = 1
  ) +
  scale_fill_manual(values = pred_cols) +
  scale_color_manual(values = pred_cols) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(x = "Predator Presence", y = "β (density‐dependence)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

ggsave(
  here::here("figures", "predators_vs_beta_boxplot.png"),
  plot = p1, width = 5, height = 5, dpi = 300
)


# ──────────────────────────────────────────────────────────────────────
# 8. Paired‐vs‐unpaired study plot (revised styling)
# ──────────────────────────────────────────────────────────────────────

library(ggplot2)
library(scales)   # for alpha()

# unpaired_studies and paired_studies already exist from above
# paired_cols is your named vector of colors by paired_substudy_num

# create a semi‐transparent fill vector for paired points
paired_fill <- setNames(alpha(paired_cols, 0.7), names(paired_cols))

p2 <- ggplot() +
  # 0) Zero reference line
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "darkgray",
    size       = 0.8
  ) +
  # 1) Unpaired studies: light gray bubbles
  geom_jitter(
    data   = unpaired_studies,
    aes(x = predators, y = sinh(betanls2_asinh)),
    shape  = 21, width = 0.1, size = 2.5,
    fill   = "gray85", color = "gray50",
    stroke = 0.6, alpha = 0.8
  ) +
  # 2) Paired studies: colored connecting lines
  geom_path(
    data   = paired_studies,
    aes(
      x     = predators,
      y     = sinh(betanls2_asinh),
      group = paired_substudy_num,
      color = paired_substudy_num
    ),
    size = 1
  ) +
  # 3) Paired study points: semi‐opaque fill + black outline
  geom_point(
    data   = paired_studies,
    aes(
      x    = predators,
      y    = sinh(betanls2_asinh),
      fill = paired_substudy_num
    ),
    shape  = 21, size = 2,
    color  = "black", stroke = 0.8
  ) +
  # 4) Meta‐prediction CIs: thick black bars with caps
  geom_errorbar(
    data    = new_pred,
    aes(
      x    = predators,
      y    = beta_hat,
      ymin = ci_lo,
      ymax = ci_hi
    ),
    width   = 0.15,
    size    = 1,
    color   = "black",
    lineend = "round"
  ) +
  # 5) Meta‐prediction points: large diamond, bold outline
  geom_point(
    data   = new_pred,
    aes(x = predators, y = beta_hat),
    shape  = 23, size = 4,
    fill   = "#1E3A5F",  # or your chosen color
    color  = "black", stroke = 1.2
  ) +
  # scales for paired studies
  scale_color_manual(values = paired_cols, guide = "none") +
  scale_fill_manual(values = paired_fill, guide = "none") +
  # axes
  scale_x_discrete(labels = c("Absent", "Present")) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  # proper axis labels
  labs(
    x = "Predator Presence",
    y = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", cm^2, ~ fish^-1, ~ day^-1, ")"
      )
    )
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

print(p2)

ggsave(
  here::here("figures", "Fig4_paired_vs_unpaired.png"),
  plot   = p2,
  width  = 6, height = 6,
  dpi    = 300,
  bg     = "white"
)
