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
  
  # LOG10 x–axis with 10^x tick labels:
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  
  labs(
    x = "Mean Density",
    y = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", 
        cm^2,   # cm²
        ~ fish^-1, 
        ~ day^-1,
        ")"
      )
    )
  ) +
  theme_classic(base_size = 10)

print(p1)

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
    y     = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", 
        cm^2,   # cm²
        ~ fish^-1, 
        ~ day^-1,
        ")"
      )
    )
    # title = "Duration vs. β"
  ) +
  theme_classic(base_size = 10)

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
    y     = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", 
        cm^2,   # cm²
        ~ fish^-1, 
        ~ day^-1,
        ")"
      )
    )
    ) +
  theme_classic(base_size = 10)

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
    y     = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", 
        cm^2,   # cm²
        ~ fish^-1, 
        ~ day^-1,
        ")"
      )
    ) ) +
  theme_classic(base_size = 10)

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
    y     = expression(
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
    title = "Study Type vs. β"
  ) +
  theme_classic(base_size = 10) +
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
  
  # axis labels using expression()
  labs(
    x = expression(paste(
      "Body Size at Experiment Start (", mm, ")"
    )),
    y = expression(paste(
      "Strength of density-dependent mortality, ",
      beta,
      " (", cm^2,  ~ fish^-1, ~ day^-1, ")"
    ))
  ) +
  
  theme_classic(base_size = 10) +
  theme(
    strip.text       = element_text(face = "italic", size = 10),
    axis.title       = element_text(size = 16, face = "bold"),
    axis.text        = element_text(size = 14),
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
# 10. Experiment Duration by Species (selected panels)
# ----------------------------------------------------------------------------

# (a) List of g_sp to plot
sel_sp <- c(
  "Pomacentrus_amboinensis",
  "Dascyllus_aruanus",
  "Coryphopterus_glaucofraenum",
  "Rhinogobiops_nicholsii",
  "Thalassoma_bifasciatum",
  "Stegastes_partitus",
  "Elacatinus_evelynae",
  "Halichoeres_garnoti",
  "Chrysiptera_parasema",
  "Trachinops_caudimaculatus",
  "Dascyllus_melanurus"
)

# (b) Subset and prepare labels
dur_sp_df <- all_dat2 %>%
  filter(
    predators == "present",
    g_sp %in% sel_sp
  ) %>%
  mutate(
    beta_raw = betanls2_raw_cm,
    Species  = gsub("_", " ", g_sp)  # convert to normal spacing
  )

# (c) Make the 11‐panel duration vs beta plot
p10 <- ggplot(dur_sp_df, aes(x = duration, y = beta_raw)) +
  geom_point(
    color = "#4682B4",
    size  = 2.5,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se     = FALSE,
    color  = "black",
    size   = 1
  ) +
  facet_wrap(
    ~ Species,
    scales = "free",
    ncol   = 4
  ) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(
    x = "Experiment Duration (Days)",
    y = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", cm^2,  ~ fish^-1,  ~ day^-1,  ")"
      )
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text        = element_text(face = "italic", size = 8),
    axis.title        = element_text(size = 10, face = "bold"),
    axis.text         = element_text(size = 10),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(5, 5, 5, 5)
  )

print(p10)

# (d) Save
ggsave(
  filename = here::here("figures", "duration_by_species.png"),
  plot     = p10,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)


