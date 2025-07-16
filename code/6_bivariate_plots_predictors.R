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


# (a) Helper: the factor levels for expt_obs (so “Exp” is baseline)
exp_levels <- levels(df$expt_obs)

# (b) Density curve over mean_density
logmd_seq <- seq(min(df$logmd_c), max(df$logmd_c), length.out = 100)
df_den <- tibble(
  logmd_c  = logmd_seq,
  dur_c    = 0,
  size_c   = 0,
  max_c    = 0,
  expt_obs = factor("Exp", levels = exp_levels)
)
X_den      <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_den)
newmods_den <- X_den[, -1]  # drop intercept
df_den <- df_den %>%
  bind_cols(predict_beta(m_all_scaled, newmods_den)) %>%
  mutate(mean_density = exp(logmd_c * sd_logden + mean_logden))

# (c) Duration curve
dur_seq <- seq(min(df$dur_c), max(df$dur_c), length.out = 100)
df_dur <- tibble(
  logmd_c  = 0,
  dur_c    = dur_seq,
  size_c   = 0,
  max_c    = 0,
  expt_obs = factor("Exp", levels = exp_levels)
)
X_dur       <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_dur)
newmods_dur <- X_dur[, -1]
df_dur <- df_dur %>%
  bind_cols(predict_beta(m_all_scaled, newmods_dur)) %>%
  mutate(duration = dur_c * sd_dur + mean_dur)

# (d) Maximum length curve
max_seq <- seq(min(df$max_c), max(df$max_c), length.out = 100)
df_max <- tibble(
  logmd_c  = 0,
  dur_c    = 0,
  size_c   = 0,
  max_c    = max_seq,
  expt_obs = factor("Exp", levels = exp_levels)
)
X_max       <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_max)
newmods_max <- X_max[, -1]
df_max <- df_max %>%
  bind_cols(predict_beta(m_all_scaled, newmods_max)) %>%
  mutate(max_length_cm = max_c * sd_max + mean_max)

# (e) Initial size curve
size_seq <- seq(min(df$size_c), max(df$size_c), length.out = 100)
df_size <- tibble(
  logmd_c  = 0,
  dur_c    = 0,
  size_c   = size_seq,
  max_c    = 0,
  expt_obs = factor("Exp", levels = exp_levels)
)
X_size        <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_size)
newmods_size  <- X_size[, -1]
df_size <- df_size %>%
  bind_cols(predict_beta(m_all_scaled, newmods_size)) %>%
  mutate(size_start = size_c * sd_size + mean_size)

# (f) Experimental vs Observational predictions
df_eo <- tibble(
  logmd_c  = 0,
  dur_c    = 0,
  size_c   = 0,
  max_c    = 0,
  expt_obs = factor(c("Exp","Obs"), levels = exp_levels)
)
X_eo       <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_eo)
newmods_eo <- X_eo[, -1]
df_eo <- df_eo %>% bind_cols(predict_beta(m_all_scaled, newmods_eo))


# --- now your existing p1–p5 code will find df_den, df_dur, df_max, df_size, df_eo ---

## (1) Density
p1 <- ggplot(df_den, aes(mean_density, beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4682B4", alpha = 0.2) +
  geom_line(color = "#4682B4", size = 1) +
  geom_point(data = df, aes(mean_density, betanls2_raw_cm),
             color = "#4682B4", alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    name   = expression(Mean~density~(fish~m^{-2}))
  ) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000,-100,-10,-1,0,1,10,100,1000),
    name   = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", cm^2, ~ fish^-1, ~ day^-1, ")"
      )
    )
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6, color = "black"),
    plot.margin = margin(2,2,2,2)
  )

## (2) Duration
p2 <- ggplot(df_dur, aes(duration, beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4682B4", alpha = 0.2) +
  geom_line(color = "#4682B4", size = 1) +
  geom_point(data = df, aes(duration, betanls2_raw_cm),
             color = "#4682B4", alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Duration (days)", y = NULL) +
  scale_y_continuous(trans = "asinh",
                     breaks = c(-1000,-100,-10,-1,0,1,10,100,1000),
                     name = expression(
                       paste(
                         "Strength of density-dependent mortality, ",
                         beta,
                         " (", cm^2, ~ fish^-1, ~ day^-1, ")"
                       )
                     )) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.x = element_text(size = 8),
    axis.text    = element_text(size = 6, color = "black"),
    plot.margin  = margin(2,2,2,2)
  )

## (3) Max Length
p3 <- ggplot(df_max, aes(max_length_cm, beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4682B4", alpha = 0.2) +
  geom_line(color = "#4682B4", size = 1) +
  geom_point(data = df, aes(max_length_cm, betanls2_raw_cm),
             color = "#4682B4", alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Maximum Length (cm)", y = NULL) +
  scale_y_continuous(trans = "asinh",
                     breaks = c(-1000,-100,-10,-1,0,1,10,100,1000),
                     name = expression(
                       paste(
                         "Strength of density-dependent mortality, ",
                         beta,
                         " (", cm^2, ~ fish^-1, ~ day^-1, ")"
                       )
                     )) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.x = element_text(size = 8),
    axis.text    = element_text(size = 6, color = "black"),
    plot.margin  = margin(2,2,2,2)
  )

## (4) Initial Size
p4 <- ggplot(df_size, aes(size_start, beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4682B4", alpha = 0.2) +
  geom_line(color = "#4682B4", size = 1) +
  geom_point(data = df, aes(size_start, betanls2_raw_cm),
             color = "#4682B4", alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Initial Size (cm)", y = NULL) +
  scale_y_continuous(trans = "asinh",
                     breaks = c(-1000,-100,-10,-1,0,1,10,100,1000),
                     name = expression(
                       paste(
                         "Strength of density-dependent mortality, ",
                         beta,
                         " (", cm^2, ~ fish^-1, ~ day^-1, ")"
                       )
                     )) +
  theme_classic(base_size = 10) +
  theme(
    axis.title.x = element_text(size = 10),
    axis.text    = element_text(size = 10, color = "black"),
    plot.margin  = margin(2,2,2,2)
  )

## (5) Study Type
p5 <- ggplot() +
  geom_jitter(
    data  = df,
    aes(x = expt_obs, y = betanls2_raw_cm, fill = expt_obs),
    shape = 21, color = "black", width = 0.15, size = 1.5, alpha = 0.6
  ) +
  geom_errorbar(
    data  = df_eo,
    aes(x = expt_obs, ymin = lower, ymax = upper),
    width = 0.1, size = 0.8, color = "black"
  ) +
  geom_point(
    data = df_eo,
    aes(x = expt_obs, y = beta),
    shape = 21, fill = "black", size = 3
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Study Type", y = expression(
    paste(
      "Strength of density-dependent mortality, ",
      beta,
      " (", cm^2, ~ fish^-1, ~ day^-1, ")"
    )
  )) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000,-100,-10,-1,0,1,10,100,1000)
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6, color = "black"),
    legend.position = "none",
    plot.margin     = margin(2,2,2,2)
  )

# ----------------------------------------------------------------------------
# 7. Display & Save
# ----------------------------------------------------------------------------
print(p1); print(p2); print(p3); print(p4); print(p5)

ggsave(here("figures","bivar_density_vs_beta.png"),      p1, width=5, height=4, dpi=300)
ggsave(here("figures","bivar_duration_vs_beta.png"),     p2, width=5, height=4, dpi=300)
ggsave(here("figures","bivar_maxlen_vs_beta.png"),       p3, width=5, height=4, dpi=300)
ggsave(here("figures","bivar_sizestart_vs_beta.png"),    p4, width=5, height=4, dpi=300)
ggsave(here("figures","bivar_exptobs_vs_beta.png"),      p5, width=5, height=4, dpi=300)
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



# ----------------------------------------------------------------------------
# 8. Body Size by Species
# ----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

# 8a. Filter and prepare data
all_filtered <- all_dat2 %>%
  filter(predators == "present") %>%
  select(
    size_start, 
    betanls2_raw_cm, 
    study_num, 
    substudy_num, 
    g_sp, 
    family
  ) %>%
  filter(!substudy_num %in% c(119, 121)) %>%
  drop_na() %>%
  group_by(g_sp) %>%
  filter(n() > 3,
         !g_sp %in% c(
           "Dascyllus_flavicaudus",
           "Sebastes_atrovirens",
           "Thalassoma_hardwicke"
         )
  ) %>%
  ungroup() %>%
  mutate(
    # convert underscores to spaces and preserve order
    g_sp = str_replace_all(g_sp, "_", " "),
    g_sp = factor(g_sp, levels = unique(g_sp))
  )

# 8b. Plot: body size vs. β, facetted by species
p_body_size <- ggplot(all_filtered, aes(x = size_start, y = betanls2_raw_cm)) +
  geom_point(
    color = "black", fill = "steelblue", shape = 21, 
    alpha = 0.6, size = 3
  ) +
  stat_smooth(
    method = "lm", formula = y ~ x, 
    se = FALSE, color = "black", size = 1
  ) +
  facet_wrap(~ g_sp, scales = "free", ncol = 4) +
  
  # match the other panels' y‐axis formatting
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
    labels = comma
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "gray"
  ) +
  
  # axis labels using expression()
  labs(
    x = expression(paste(
      "Initial Size (", mm, ")"
    )),
    y = expression(paste(
      "Strength of density-dependent mortality, ",
      beta,
      " (", cm^2,  ~ fish^-1, ~ day^-1, ")"
    ))
  ) +
  
  theme_classic(base_size = 10) +
  theme(
    strip.text       = element_text(face = "italic", size = 8),
    axis.title       = element_text(size = 12, face = "bold"),
    axis.text        = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

ggsave(
  filename = here::here("figures", "body_size_by_species.png"),
  plot     = p_body_size,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)



# ----------------------------------------------------------------------------
# 8b. Experimental Duration by Species
# ----------------------------------------------------------------------------

# (a) Filter and prepare data — same species set as 8a, but include duration
dur_filtered <- all_dat2 %>%
  filter(predators == "present") %>%
  select(
    duration,
    betanls2_raw_cm,
    study_num,
    substudy_num,
    g_sp
  ) %>%
  filter(!substudy_num %in% c(119, 121)) %>%
  drop_na() %>%
  group_by(g_sp) %>%
  filter(
    n() > 3,
    !g_sp %in% c(
      "Dascyllus_flavicaudus",
      "Sebastes_atrovirens",
      "Thalassoma_hardwicke"
    )
  ) %>%
  ungroup() %>%
  mutate(
    # convert underscores to spaces and preserve order
    g_sp = str_replace_all(g_sp, "_", " "),
    g_sp = factor(g_sp, levels = unique(g_sp))
  )

# (b) Plot: duration vs. β, facetted by species
p_duration_by_species <- ggplot(dur_filtered, aes(x = duration, y = betanls2_raw_cm)) +
  geom_point(
    color = "black", fill = "steelblue", shape = 21,
    alpha = 0.6, size = 3
  ) +
  stat_smooth(
    method = "lm", formula = y ~ x,
    se = FALSE, color = "black", size = 1
  ) +
  facet_wrap(~ g_sp, scales = "free", ncol = 4) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
    labels = comma
  ) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "gray"
  ) +
  labs(
    x = "Experiment Duration (Days)",
    y = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", cm^2,  ~ fish^-1, ~ day^-1, ")"
      )
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text       = element_text(face = "italic", size = 8),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

# (c) Save
ggsave(
  filename = here("figures", "duration_by_species.png"),
  plot     = p_duration_by_species,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)


# ----------------------------------------------------------------------------
# 9. Experimental vs. Observational: Paired & Unpaired Studies
# ----------------------------------------------------------------------------

# (a) Define palette for species (if not already defined)
pal2 <- wesanderson::wes_palette("Zissou1", 12, type = "continuous")

# (b) Subset to predator-present studies with a defined expt_obs
expobs <- all_dat2 %>%
  filter(
    predators == "present",
    !is.na(expt_obs)
  ) %>%
  mutate(
    expt_obs = factor(expt_obs, levels = c("Exp", "Obs"))
  )

# (c) Split into paired vs. unpaired
paired_expobs   <- expobs %>% filter(!is.na(expt_obs_pairs))
unpaired_expobs <- expobs %>% filter(is.na(expt_obs_pairs))

# (d) Build meta-analytic means & CIs from the full model predictions (df_eo)
eo_coef <- df_eo %>%
  transmute(
    expobs   = expt_obs,
    betanls2 = beta,
    lower    = lower,
    upper    = upper
  )

# (e) Plot paired and unpaired studies with meta-analytic estimates
p_expobs <- ggplot() +
  # Unpaired data as jittered points
  geom_jitter(
    data    = unpaired_expobs,
    aes(x = expt_obs, y = sinh(betanls2_asinh)),
    color   = "black", alpha = 0.5, width = 0.1
  ) +
  # Paired data with connecting paths
  geom_path(
    data    = paired_expobs,
    aes(
      x     = expt_obs,
      y     = sinh(betanls2_asinh),
      group = expt_obs_pairs,
      color = g_sp
    ),
    linewidth = 0.8, alpha = 0.6
  ) +
  geom_point(
    data  = paired_expobs,
    aes(
      x     = expt_obs,
      y     = sinh(betanls2_asinh),
      color = g_sp
    ),
    size  = 3, alpha = 0.8
  ) +
  # Meta-analytic means and CIs
  geom_errorbar(
    data    = eo_coef,
    aes(
      x    = expobs,
      y    = betanls2,
      ymin = lower,
      ymax = upper
    ),
    width     = 0.2,
    color     = "red",
    linewidth = 0.8
  ) +
  geom_point(
    data  = eo_coef,
    aes(x = expobs, y = betanls2),
    color = "red",
    shape = 18,
    size  = 4
  ) +
  # Zero reference line
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "gray"
  ) +
  # Sinh-transformed y-axis
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  scale_color_manual(
    values = pal2,
    name   = "Species"
  ) +
  labs(
    x = "Study Type",
    y = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", cm^2, ~ fish^-1, ~ day^-1, ")"
      )
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text      = element_text(size = 10),
    axis.title     = element_text(size = 10),
    legend.position = "right"
  )

print(p_expobs)

# (f) Save to file
ggsave(
  filename = here::here("figures", "exp_obs_paired_vs_unpaired.png"),
  plot     = p_expobs,
  width    = 6,
  height   = 4,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)



# ----------------------------------------------------------------------------
# 10. Combined bivariate panels: Density, Duration & Max Length
# ----------------------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(here)

# 10a) Shared components ----------------------------------------------------
hline0 <- geom_hline(yintercept = 0, linetype = "dashed", color = "black")

y_scale <- scale_y_continuous(
  trans  = "asinh",
  breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
  name   = expression(
    paste(
      "Strength of density-dependent mortality, ",
      beta,
      " (", cm^2, ~ fish^-1, ~ day^-1, ")"
    )
  )
)

theme_bivar <- theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text    = element_text(size = 10, color = "black"),
    plot.margin  = margin(3, 3, 3, 3)
  )

# 10b) Panel A: Mean Density -----------------------------------------------
p1 <- ggplot() +
  geom_ribbon(
    data = df_den,
    aes(x = mean_density, y = beta, ymin = lower, ymax = upper),
    fill  = "#4682B4", alpha = 0.2
  ) +
  geom_line(
    data = df_den,
    aes(x = mean_density, y = beta),
    color = "#4682B4", size = 1
  ) +
  geom_point(
    data = df,
    aes(x = mean_density, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  hline0 +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    name   = expression(Mean~density~(fish~m^{-2}))
  ) +
  y_scale +
  theme_bivar

# 10c) Panel B: Duration ----------------------------------------------------
p2 <- ggplot() +
  geom_ribbon(
    data = df_dur,
    aes(x = duration, y = beta, ymin = lower, ymax = upper),
    fill  = "#4682B4", alpha = 0.2
  ) +
  geom_line(
    data = df_dur,
    aes(x = duration, y = beta),
    color = "#4682B4", size = 1
  ) +
  geom_point(
    data = df,
    aes(x = duration, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  hline0 +
  labs(x = "Duration (days)") +
  y_scale +
  theme_bivar +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# 10d) Panel C: Max Length -------------------------------------------------
p3 <- ggplot() +
  geom_ribbon(
    data = df_max,
    aes(x = max_length_cm, y = beta, ymin = lower, ymax = upper),
    fill  = "#4682B4", alpha = 0.2
  ) +
  geom_line(
    data = df_max,
    aes(x = max_length_cm, y = beta),
    color = "#4682B4", size = 1
  ) +
  geom_point(
    data = df,
    aes(x = max_length_cm, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  hline0 +
  labs(x = "Maximum Length (cm)") +
  y_scale +
  theme_bivar +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# ----------------------------------------------------------------------------
# 10. Combined bivariate panels: Density, Duration & Max Length (cowplot)
# ----------------------------------------------------------------------------
library(ggplot2)
library(cowplot)
library(here)

# Assume p1, p2, p3 are defined as in the previous steps.

# Add a consistent margin so panel content doesn’t shift
p1m <- p1 + theme(plot.margin = margin(5,5,5,5))
p2m <- p2 + theme(plot.margin = margin(5,5,5,5))
p3m <- p3 + theme(plot.margin = margin(5,5,5,5))

# Combine and align
combined <- plot_grid(
  p1m, p2m, p3m,
  nrow           = 1,
  align          = "v",   # align vertically
  axis           = "l",   # align on left axes
  rel_widths     = c(1,1,1),
  labels         = c("A", "B", "C"),
  label_size     = 14,
  label_fontface = "bold",
  # place labels at 90% from left, 95% from bottom
  label_x        = c(0.90, 0.90, 0.90),
  label_y        = c(0.95, 0.95, 0.95)
)

# Render & Save
print(combined)

ggsave(
  filename = here("figures", "bivar_density_duration_maxlen.pdf"),
  plot     = combined,
  width    = 12,
  height   = 4,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)
