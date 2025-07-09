################################################################################
# File: 6_bivariate_plots.R
# Purpose: Generate bivariate plots of each predictor vs. β, holding all other
#          covariates at their mean (i.e., zero on the centered scale)
# Author: Adrian Stier (rev. 2025-07-09)
################################################################################

# ----------------------------------------------------------------------------
# 1. Dependencies & Data Load
# ----------------------------------------------------------------------------
source(here::here("code", "0_libraries.R"))               # Core packages
source(here::here("code", "1_data_phylogeny_loading.R"))  # Loads all_dat2, pruned_tree, etc.

# ensure output directory exists
dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)

library(dplyr)    # data manipulation
library(ggplot2)  # plotting
library(metafor)  # meta-analysis
library(tibble)   # tibble utilities

# ----------------------------------------------------------------------------
# 2. Data Preparation (centered covariates)
# ----------------------------------------------------------------------------
df <- all_dat2 %>%
  filter(predators == "present") %>%
  mutate(
    expt_obs = factor(expt_obs, levels = c("Exp", "Obs")),
    logmd_c  = as.numeric(scale(log(mean_density))),
    dur_c    = as.numeric(scale(duration)),
    size_c   = as.numeric(scale(size_start)),
    max_c    = as.numeric(scale(max_length_cm))
  )

phylo_vcv <- ape::vcv(pruned_tree, corr = TRUE)
keep_sp    <- intersect(rownames(phylo_vcv), unique(df$g_sp))
phylo_vcv  <- phylo_vcv[keep_sp, keep_sp]
df         <- df %>% filter(g_sp %in% keep_sp)

rand_list  <- list(~1 | study_num/substudy_num, ~1 | g_sp)

# ----------------------------------------------------------------------------
# 3. Helper: predict & back-transform β
# ----------------------------------------------------------------------------
predict_beta <- function(model, newmods) {
  pr <- predict(model, newmods = newmods)
  tibble(
    beta  = sinh(pr$pred),
    lower = sinh(pr$ci.lb),
    upper = sinh(pr$ci.ub)
  )
}

# ----------------------------------------------------------------------------
# 4. Fit the no-interaction additive model
# ----------------------------------------------------------------------------
m_all_scaled <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ logmd_c + dur_c + size_c + max_c + expt_obs,
  random = rand_list,
  R      = list(g_sp = phylo_vcv),
  data   = df,
  method = "REML",
  test   = "t"
)

# ----------------------------------------------------------------------------
# 5. Compute means & SDs for un-centering
# ----------------------------------------------------------------------------
mean_logden <- mean(log(df$mean_density), na.rm = TRUE)
sd_logden   <- sd(  log(df$mean_density), na.rm = TRUE)
mean_dur    <- mean(df$duration,      na.rm = TRUE)
sd_dur      <- sd(  df$duration,      na.rm = TRUE)
mean_size   <- mean(df$size_start,    na.rm = TRUE)
sd_size     <- sd(  df$size_start,    na.rm = TRUE)
mean_max    <- mean(df$max_length_cm, na.rm = TRUE)
sd_max      <- sd(  df$max_length_cm, na.rm = TRUE)

# ----------------------------------------------------------------------------
# 6. Bivariate curves & raw data points
# ----------------------------------------------------------------------------

## (1) Density (log scale)
logmd_c_seq <- seq(min(df$logmd_c), max(df$logmd_c), length.out = 120)
new_den <- tibble(
  logmd_c  = logmd_c_seq,
  dur_c    = 0,
  size_c   = 0,
  max_c    = 0,
  expt_obs = factor("Exp", levels = levels(df$expt_obs))
)
X_den   <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, data = new_den)[,-1]
pred_den <- predict_beta(m_all_scaled, X_den)
df_den   <- bind_cols(new_den, pred_den) %>%
  mutate(mean_density = exp(logmd_c * sd_logden + mean_logden))

p1 <- ggplot(df_den, aes(x = mean_density, y = beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1.2) +
  geom_point(
    data = df,
    aes(x = mean_density, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(
    trans = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(
    x     = "Mean Density (log scale)",
    y     = expression(beta ~ ": density-dependent mortality"),
    title = "Density vs. β"
  ) +
  theme_classic(base_size = 14)

## (2) Duration
dur_c_seq <- seq(min(df$dur_c), max(df$dur_c), length.out = 120)
new_dur   <- tibble(
  logmd_c  = 0,
  dur_c    = dur_c_seq,
  size_c   = 0,
  max_c    = 0,
  expt_obs = factor("Exp", levels = levels(df$expt_obs))
)
X_dur    <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, data = new_dur)[,-1]
pred_dur <- predict_beta(m_all_scaled, X_dur)
df_dur   <- bind_cols(new_dur, pred_dur) %>%
  mutate(duration = dur_c * sd_dur + mean_dur)

p2 <- ggplot(df_dur, aes(x = duration, y = beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1.2) +
  geom_point(
    data = df,
    aes(x = duration, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(
    trans = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(
    x     = "Duration (days)",
    y     = expression(beta ~ ": density-dependent mortality"),
    title = "Duration vs. β"
  ) +
  theme_classic(base_size = 14)

## (3) Max Length
max_c_seq <- seq(min(df$max_c), max(df$max_c), length.out = 120)
new_max   <- tibble(
  logmd_c  = 0,
  dur_c    = 0,
  size_c   = 0,
  max_c    = max_c_seq,
  expt_obs = factor("Exp", levels = levels(df$expt_obs))
)
X_max    <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, data = new_max)[,-1]
pred_max <- predict_beta(m_all_scaled, X_max)
df_max   <- bind_cols(new_max, pred_max) %>%
  mutate(max_length_cm = max_c * sd_max + mean_max)

p3 <- ggplot(df_max, aes(x = max_length_cm, y = beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1.2) +
  geom_point(
    data = df,
    aes(x = max_length_cm, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(
    trans = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(
    x     = "Max Length (cm)",
    y     = expression(beta ~ ": density-dependent mortality"),
    title = "Max Length vs. β"
  ) +
  theme_classic(base_size = 14)

## (4) Initial Size
size_c_seq <- seq(min(df$size_c), max(df$size_c), length.out = 120)
new_size   <- tibble(
  logmd_c  = 0,
  dur_c    = 0,
  size_c   = size_c_seq,
  max_c    = 0,
  expt_obs = factor("Exp", levels = levels(df$expt_obs))
)
X_size    <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, data = new_size)[,-1]
pred_size <- predict_beta(m_all_scaled, X_size)
df_size   <- bind_cols(new_size, pred_size) %>%
  mutate(size_start = size_c * sd_size + mean_size)

p4 <- ggplot(df_size, aes(x = size_start, y = beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1.2) +
  geom_point(
    data = df,
    aes(x = size_start, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(
    trans = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(
    x     = "Initial Size (cm)",
    y     = expression(beta ~ ": density-dependent mortality"),
    title = "Initial Size vs. β"
  ) +
  theme_classic(base_size = 14)

## (5) Study Type (Exp vs Obs)
new_eo   <- tibble(
  logmd_c  = 0,
  dur_c    = 0,
  size_c   = 0,
  max_c    = 0,
  expt_obs = factor(c("Exp", "Obs"), levels = levels(df$expt_obs))
)
X_eo     <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, data = new_eo)[,-1]
pred_eo  <- predict_beta(m_all_scaled, X_eo)
df_eo    <- bind_cols(new_eo, pred_eo)

p5 <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = expt_obs, y = betanls2_raw_cm, fill = expt_obs),
    shape = 21, color = "black", width = 0.2, size = 2, alpha = 0.6
  ) +
  geom_errorbar(
    data = df_eo,
    aes(x = expt_obs, ymin = lower, ymax = upper),
    width = 0.1, size = 1.2
  ) +
  geom_point(
    data = df_eo,
    aes(x = expt_obs, y = beta),
    shape = 21, fill = "black", size = 4
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(
    trans = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(
    x     = "Study Type",
    y     = expression(beta ~ ": density-dependent mortality"),
    title = "Study Type vs. β"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# ----------------------------------------------------------------------------
# 7. Display & Save
# ----------------------------------------------------------------------------
print(p1); print(p2); print(p3); print(p4); print(p5)

ggsave(
  here::here("figures", "bivar_density_vs_beta.png"),
  p1, width = 5, height = 4, dpi = 300
)
ggsave(
  here::here("figures", "bivar_duration_vs_beta.png"),
  p2, width = 5, height = 4, dpi = 300
)
ggsave(
  here::here("figures", "bivar_maxlen_vs_beta.png"),
  p3, width = 5, height = 4, dpi = 300
)
ggsave(
  here::here("figures", "bivar_sizestart_vs_beta.png"),
  p4, width = 5, height = 4, dpi = 300
)
ggsave(
  here::here("figures", "bivar_exptobs_vs_beta.png"),
  p5, width = 5, height = 4, dpi = 300
)