###############################################################################
# File: 5_model_selection_and_autocorrelation.R
# Purpose: 
#   1. Explore autocorrelation among predictors  
#   2. Visualize β vs. density by study type  
#   3. Prepare data & phylogeny  
#   4. Fit full 5-way model and two-way interaction model  
#   5. Stepwise removal by AICc and by p-value  
#   6. Summarize best-fitting additive model  
# Author: Adrian Stier (rev. 2025-07-09)
###############################################################################

# ----------------------------------------------------------------------------
# 1. Dependencies & Data Load
# ----------------------------------------------------------------------------
source(here::here("code", "0_libraries.R"))               # Core packages
source(here::here("code", "1_data_phylogeny_loading.R"))  # Loads all_dat2, pruned_tree, etc.

# ensure output directories exist
dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("results"), showWarnings = FALSE, recursive = TRUE)

library(GGally)    # for ggpairs
library(ggplot2)   # plotting
library(metafor)   # meta-analysis
library(dplyr)     # data manipulation
library(gt)        # tables
library(tibble)    # rownames_to_column
library(broom)     # tidy rma.mv objects

# ----------------------------------------------------------------------------
# 2. Autocorrelation among predictor variables (exploratory)
# ----------------------------------------------------------------------------
all <- all_dat2  # temporary alias

# compute log density
all$logmeandensity <- log(all$mean_density)

# select raw predictors
predictors <- all %>%
  select(expt_obs, size_start, duration, logmeandensity, max_length_cm, max_weight_g)

# multipanel scatterplot matrix
ggpairs(
  predictors,
  lower = list(continuous = wrap("smooth", color = "blue")),
  diag  = list(continuous = "densityDiag"),
  upper = list(continuous = wrap("cor", size = 5))
) +
  theme_minimal()

# ----------------------------------------------------------------------------
# 3. β vs. density by study type (exploratory)
# ----------------------------------------------------------------------------
log_ticks  <- c(1e-4,1e-3,1e-2,1e-1,1,10,100,1000)
sinh_ticks <- sinh(c(-10,-5,-2,-1,0,1,2,5,10))

ggplot(all, aes(mean_density, betanls2_raw_cm, shape=expt_obs, color=expt_obs)) +
  geom_point(alpha=0.5, size=2) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(trans="log10", breaks=log_ticks, labels=log_ticks) +
  scale_y_continuous(trans="asinh", breaks=sinh_ticks) +
  scale_shape_manual(values=c(16,17)) +
  scale_color_manual(values=c("blue","red")) +
  labs(
    x     = "Mean Density (n₀/m², log scale)",
    y     = "Density-dependence strength (β)",
    title = "β vs. Density by Study Type"
  ) +
  theme_classic(base_size=14) +
  theme(legend.position="top")

# ----------------------------------------------------------------------------
# 4. Data preparation & phylogeny
# ----------------------------------------------------------------------------
all <- all_dat2 %>%
  filter(predators == "present") %>%
  mutate(
    expt_obs = factor(expt_obs, levels=c("Exp","Obs")),
    logmd_c  = as.numeric(scale(log(mean_density))),
    dur_c    = as.numeric(scale(duration)),
    size_c   = as.numeric(scale(size_start)),
    max_c    = as.numeric(scale(max_length_cm))
  )

phylo_vcv <- ape::vcv(pruned_tree, corr=TRUE)
keep_sp    <- intersect(rownames(phylo_vcv), unique(all$g_sp))
phylo_vcv  <- phylo_vcv[keep_sp, keep_sp]
all         <- all %>% filter(g_sp %in% keep_sp)

rand_list <- list(~1 | study_num/substudy_num, ~1 | g_sp)

# ----------------------------------------------------------------------------
# 5. Helper functions: fit rma.mv & AICc
# ----------------------------------------------------------------------------
fit_rma <- function(formula_obj) {
  rma.mv(
    yi     = betanls2_asinh,
    V      = betanlsvar_asinh,
    mods   = formula_obj,
    random = rand_list,
    R      = list(g_sp = phylo_vcv),
    data   = all,
    method = "REML",
    test   = "t"
  )
}

get_AICc <- function(mod) {
  summary(mod)$fit.stats["AICc","REML"] %>% as.numeric()
}

# ----------------------------------------------------------------------------
# 6. Build “full” 5-way model & record AICc
# ----------------------------------------------------------------------------
vars_centered <- c("expt_obs","logmd_c","dur_c","size_c","max_c")
all_terms     <- unlist(lapply(seq_along(vars_centered), function(k)
  combn(vars_centered, k, FUN=paste, collapse=":", simplify=TRUE)))
full_formula  <- reformulate(all_terms)

m_current     <- fit_rma(full_formula)
current_AICc  <- get_AICc(m_current)
message("Initial AICc (full model) = ", round(current_AICc,3))

# ----------------------------------------------------------------------------
# 7. Support functions for hierarchical deletion
# ----------------------------------------------------------------------------
is_nested_in_higher <- function(lower, terms) {
  parts <- strsplit(lower,":",fixed=TRUE)[[1]]
  any(sapply(terms, function(cand) {
    cp <- strsplit(cand,":",fixed=TRUE)[[1]]
    length(cp)>length(parts) && all(parts %in% cp)
  }))
}

get_coef_rows <- function(label, names) {
  parts <- strsplit(label,":",fixed=TRUE)[[1]]
  if (length(parts)==1 && parts[1]=="expt_obs")
    grep("^expt_obs(?:Exp|Obs)$", names)
  else if (parts[1]=="expt_obs") {
    sfx <- paste(parts[-1],collapse=":")
    grep(paste0("^expt_obs(?:Exp|Obs):",sfx,"$"), names)
  } else grep(paste0("^",label,"$"), names)
}

try_drop_one <- function(label, model) {
  terms_now <- attr(terms(formula(model)),"term.labels")
  parts     <- strsplit(label,":",fixed=TRUE)[[1]]
  if (length(parts)==1 || (label %in% terms_now && is_nested_in_higher(label,terms_now)))
    return(list(model=model,AICc=current_AICc,dropped=FALSE))
  rows <- get_coef_rows(label, rownames(coef(summary(model))))
  if (!length(rows)) return(list(model=model,AICc=current_AICc,dropped=FALSE))
  new_terms <- setdiff(terms_now,label)
  new_mod   <- tryCatch(fit_rma(reformulate(new_terms)), error=function(e) NULL)
  if (is.null(new_mod)) return(list(model=model,AICc=current_AICc,dropped=FALSE))
  list(model=new_mod, AICc=get_AICc(new_mod), dropped=TRUE)
}

# ----------------------------------------------------------------------------
# 8. Two-way interaction model
# ----------------------------------------------------------------------------
m_2way <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ expt_obs + logmd_c + dur_c + size_c + max_c
  + expt_obs:logmd_c + expt_obs:dur_c + expt_obs:size_c + expt_obs:max_c
  + logmd_c:dur_c   + logmd_c:size_c   + logmd_c:max_c
  + dur_c:size_c    + dur_c:max_c      + size_c:max_c,
  random = rand_list,
  R      = list(g_sp = phylo_vcv),
  data   = all,
  method = "REML",
  test   = "t"
)
summary(m_2way)

# ----------------------------------------------------------------------------
# 9. Stepwise removal: interaction-order–specific p-value thresholds
# ----------------------------------------------------------------------------
drop_if_insig <- function(label, model) {
  ct        <- coef(summary(model))
  terms_now <- attr(terms(formula(model)),"term.labels")
  parts     <- strsplit(label,":",fixed=TRUE)[[1]]
  ord       <- length(parts)
  cutoff    <- if (ord==3) 0.01 else 0.05
  
  if (ord==1 || (label %in% terms_now && is_nested_in_higher(label,terms_now)))
    return(list(model=model,dropped=FALSE))
  
  pat <- if (parts[1]=="expt_obs") {
    sfx <- sub("^expt_obs:","",label)
    paste0("^expt_obs(?:Exp|Obs):",sfx,"$")
  } else paste0("^",label,"$")
  
  idx <- grep(pat, rownames(ct))
  if (!length(idx)) return(list(model=model,dropped=FALSE))
  
  if (max(ct[idx,"pval"],na.rm=TRUE) >= cutoff) {
    new_mod <- tryCatch(fit_rma(reformulate(setdiff(terms_now,label))),
                        error=function(e) NULL)
    if (!is.null(new_mod)) return(list(model=new_mod,dropped=TRUE))
  }
  list(model=model,dropped=FALSE)
}

# [Insert hierarchical loops over 5-way→2-way terms here as in original script]

# ----------------------------------------------------------------------------
# 10. Final best-fitting additive model
# ----------------------------------------------------------------------------
m_best <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ logmd_c + dur_c + size_c + max_c + expt_obs,
  random = rand_list,
  R      = list(g_sp = phylo_vcv),
  data   = all,
  method = "REML",
  test   = "t"
)
summary(m_best)

# ----------------------------------------------------------------------------
# 11. Back‐transform & tabulate best‐fitting model
# ----------------------------------------------------------------------------
results <- broom::tidy(m_best, conf.int = TRUE) %>%
  mutate(
    estimate_bt  = sinh(estimate),
    conf.low_bt  = sinh(conf.low),
    conf.high_bt = sinh(conf.high),
    p.value      = if_else(p.value < .001, "<0.001", as.character(round(p.value, 3)))
  ) %>%
  select(
    Predictor    = term,
    Estimate     = estimate_bt,
    SE           = std.error,
    `t-value`    = statistic,
    `p-value`    = p.value,
    `CI low`     = conf.low_bt,
    `CI high`    = conf.high_bt
  )

results %>%
  gt() %>%
  tab_header(
    title    = md("**Best‐Fitting Multivariate Meta‐Analysis Model**"),
    subtitle = "Back‐Transformed Estimates"
  ) %>%
  fmt_number(columns = c("Estimate", "SE", "t-value"), decimals = 3) %>%
  cols_label(
    Predictor = "Predictor",
    SE        = "SE",
    `t-value` = "t‐value",
    `p-value` = "p‐value",
    `CI low`  = "95% CI (low)",
    `CI high` = "95% CI (high)"
  ) %>%
  tab_options(table.font.size = px(14), heading.align = "center") %>%
  gtsave(here::here("results", "best_model_summary.png"))