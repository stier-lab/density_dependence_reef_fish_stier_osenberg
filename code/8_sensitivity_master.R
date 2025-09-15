# ==============================================================================
# File: 8_sensitivity_core_beta.R
# Title: Core β Meta-Analysis Diagnostics (No Covariates)
# Author: Stier Lab (Adrian C. Stier), with assistant support
# Date: 2025-09-15
#
# PURPOSE
#   Focus exclusively on the core intercept-only meta-analysis of β, with:
#     • Heterogeneity (tau^2, sigma^2, QE/QEp, pseudo-I^2)
#     • Leave-one-out (CSV + plot)
#     • Influence diagnostics (Cook's D, hat, max|DFBETAS| + plot)
#     • Small-study effects: funnel plot + Egger test (trim-and-fill saved for SI)
#
# HOW TO RUN
#   Source this AFTER scripts 2–6 in 00_run_all.R:
#     source("8_sensitivity_core_beta.R")
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(stringr); library(ggplot2); library(forcats); library(broom); library(janitor)
  library(metafor); library(clubSandwich)
  # Optional phylogeny (off by default)
  # library(ape); library(phytools)
})

# ---------------------------- Config ------------------------------------------
ROOT_FIG <- "figs/sensitivity_core_beta"
ROOT_TAB <- "tables/sensitivity_core_beta"
dir.create(ROOT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(ROOT_TAB, recursive = TRUE, showWarnings = FALSE)

# Optional phylogeny toggle (kept for completeness; OFF by default)
USE_PHYLO <- FALSE
TREE_FILE <- "actinopt_12k_treePL.tre"  # or "1.newick.tre"

# ----------------------- Robust file location ---------------------------------
locate_file <- function(cands) {
  cands <- unique(cands)
  for (p in cands) if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
  for (pref in c(".", "output", "data", "inputs")) {
    for (p0 in cands) {
      p <- file.path(pref, p0)
      if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  if (requireNamespace("here", quietly = TRUE)) {
    for (p0 in cands) {
      p <- here::here(p0); if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
      for (pref in c("output","data","inputs")) {
        pp <- here::here(pref, p0); if (file.exists(pp)) return(normalizePath(pp, winslash = "/", mustWork = TRUE))
      }
    }
  }
  stop("Could not locate any of: ", paste(cands, collapse = ", "))
}

safe_read_csv <- function(path_or_vec) {
  if (length(path_or_vec) > 1 || !file.exists(path_or_vec)) {
    found <- locate_file(path_or_vec); message("Reading: ", found)
    readr::read_csv(found, show_col_types = FALSE)
  } else { message("Reading: ", path_or_vec); readr::read_csv(path_or_vec, show_col_types = FALSE) }
}

mk_dirs <- function(label) {
  fig_dir <- file.path(ROOT_FIG, label); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  tab_dir <- file.path(ROOT_TAB, label); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  list(fig = fig_dir, tab = tab_dir)
}
save_gg <- function(p, file, w=7, h=5, dpi=300) ggplot2::ggsave(file, p, width=w, height=h, dpi=dpi)

# --------------------------- ID + variance helpers ----------------------------
normalize_id <- function(v) {
  if (is.numeric(v)) {
    out <- ifelse(is.na(v), NA_character_, format(v, scientific = FALSE, trim = TRUE))
    sub("\\.0+$", "", out)
  } else as.character(v)
}

compute_vi_floor <- function(vi) {
  pos <- vi[is.finite(vi) & vi > 0]
  if (!length(pos)) return(NA_real_)
  max(1e-12, min(min(pos, na.rm=TRUE) * 1e-3,
                 stats::quantile(pos, probs = 0.005, na.rm = TRUE, names = FALSE)))
}

with_vi_use <- function(dat) {
  dat <- dat %>% filter(is.finite(beta), is.finite(vi))
  floor_val <- compute_vi_floor(dat$vi)
  dat %>% mutate(vi_use = if (!is.na(floor_val)) pmax(vi, floor_val) else vi,
                 vi_floor = floor_val)
}

# ----------------------------- Model helpers ----------------------------------
fit_intercept <- function(dat, method="REML", R_arg=NULL) {
  rma.mv(yi = beta, V = vi_use, mods = ~ 1,
         random = list(~ 1 | study_num/substudy_num, ~ 1 | g_sp),
         R = R_arg, data = dat, method = method, test = "t")
}

summarize_heterogeneity <- function(fit, dat) {
  tibble(
    k            = nrow(dat),
    tau2_resid   = unname(fit$tau2),
    QE           = unname(fit$QE), QEp = unname(fit$QEp),
    sigma2       = paste(round(fit$sigma2, 6), collapse = ";"),
    sigma2_names = paste(fit$s.names, collapse = ";"),
    I2_overall   = suppressWarnings(tryCatch(i2(fit)$I2, error = function(e) NA_real_))
  )
}

plot_influence <- function(infl, cooks_cutoff, fig_path) {
  dd <- tibble(
    idx     = seq_along(infl$cook.d),
    cooks_d = infl$cook.d,
    hat     = infl$hat,
    dfbetas = apply(infl$dfbs, 1, function(x) max(abs(x), na.rm = TRUE))
  )
  p <- ggplot(dd, aes(cooks_d, dfbetas)) +
    geom_point(alpha = .85) +
    geom_vline(xintercept = cooks_cutoff, linetype = "dashed") +
    labs(x = "Cook's distance", y = "max|DFBETAS|", title = "Influence diagnostics") +
    theme_minimal(base_size = 12)
  save_gg(p, fig_path, w = 7, h = 5, dpi = 300)
  dd
}

# Funnel + Egger + Trim&Fill (uses vi_use to keep consistent with fitting)
funnel_egger_trimfill <- function(dat, fig_path, egger_path, tf_path) {
  uni <- tryCatch(rma.uni(yi = beta, vi = vi_use, method = "REML", data = dat),
                  error = function(e) NULL)
  if (!is.null(uni)) {
    png(fig_path, width = 1800, height = 1400, res = 300)
    metafor::funnel(uni, main = "Funnel plot")
    dev.off()
    egger <- tryCatch(regtest(uni, model = "lm"), error=function(e) NULL)
    if (!is.null(egger)) {
      readr::write_csv(tibble(z = unname(egger$zval), p = unname(egger$pval)), egger_path)
    }
    tf <- tryCatch(trimfill(uni), error = function(e) NULL)
    if (!is.null(tf)) readr::write_csv(broom::tidy(tf), tf_path)
  }
  invisible(uni)
}

# ----------------------- Run diagnostics for baseline β -----------------------
run_core_beta <- function(dat, label = "baseline_beta", R_arg = NULL) {
  dirs <- mk_dirs(label); fig_dir <- dirs$fig; tab_dir <- dirs$tab
  
  # Per-analysis Cook's cutoff
  cooks_cutoff <- 4 / max(5, nrow(dat))
  
  # Fit
  fit <- tryCatch(fit_intercept(dat, method="REML", R_arg=R_arg), error = function(e) NULL)
  
  # Heterogeneity
  if (!is.null(fit)) {
    het <- summarize_heterogeneity(fit, dat)
    readr::write_csv(het, file.path(tab_dir, "heterogeneity.csv"))
  } else {
    readr::write_csv(tibble(note="Base intercept fit failed; heterogeneity unavailable."),
                     file.path(tab_dir, "heterogeneity.csv"))
  }
  
  # Influence + LOO (if fit exists)
  if (!is.null(fit)) {
    infl <- tryCatch(influence.rma.mv(fit), error=function(e) NULL)
    if (!is.null(infl)) {
      inf_tab <- plot_influence(infl, cooks_cutoff, file.path(fig_dir, "influence.png"))
      readr::write_csv(inf_tab %>% select(idx, cooks_d, hat, dfbetas),
                       file.path(tab_dir, "influence_metrics.csv"))
    }
    loo <- tryCatch(leave1out.rma.mv(fit), error=function(e) NULL)
    if (!is.null(loo) && length(loo$estimate) == nrow(dat)) {
      readr::write_csv(
        tibble(idx = seq_along(loo$estimate), loo_est = as.numeric(loo$estimate)),
        file.path(tab_dir, "leave_one_out.csv")
      )
      p <- ggplot(tibble(idx=seq_along(loo$estimate), y=loo$estimate),
                  aes(idx, y)) +
        geom_point(alpha=.9) +
        geom_hline(yintercept = as.numeric(fit$b), linetype="dashed") +
        labs(x="Left-out index", y="β", title="Leave-one-out β") +
        theme_minimal(12)
      save_gg(p, file.path(fig_dir, "leave_one_out.png"), w=8, h=4, dpi=300)
    }
  }
  
  # Small-study effects
  funnel_egger_trimfill(dat,
                        fig_path  = file.path(fig_dir, "funnel.png"),
                        egger_path= file.path(tab_dir, "egger_test.csv"),
                        tf_path   = file.path(tab_dir, "trimfill_SI.csv"))
  
  # Session info
  writeLines(capture.output(sessionInfo()),
             con = file.path(tab_dir, paste0("sessionInfo_", label, ".txt")))
  
  invisible(fit)
}

# =========================== Build master dataset =============================
RES_CANDS  <- c("combined_results_2024-09-17.csv",
                "output/combined_results_2024-09-17.csv")
COV_CANDS  <- c("covariates-2024-09-30.csv",
                "output/covariates-2024-09-30.csv",
                "data/covariates-2024-09-30.csv")
RAW_CANDS  <- c("all_studies_looped-2024-09-11.csv",
                "output/all_studies_looped-2024-09-11.csv",
                "data/all_studies_looped-2024-09-11.csv")

res  <- safe_read_csv(RES_CANDS) %>% clean_names() %>%
  transmute(substudy_num, beta = beta, vi = beta_variance)

covs <- safe_read_csv(COV_CANDS) %>% clean_names() %>%
  select(any_of(c("substudy_num","study_num","g_sp")))

raw  <- safe_read_csv(RAW_CANDS) %>% clean_names() %>%
  transmute(study_num, substudy_num)

# Normalize IDs to avoid join type clashes
res  <- res  %>% mutate(substudy_num = normalize_id(substudy_num))
covs <- covs %>% mutate(substudy_num = normalize_id(substudy_num),
                        study_num    = normalize_id(study_num),
                        g_sp         = as.character(g_sp))
raw  <- raw  %>% mutate(substudy_num = normalize_id(substudy_num),
                        study_num    = normalize_id(study_num))

# Merge -> analysis frame (β only)
all <- res %>%
  left_join(covs, by = "substudy_num") %>%
  left_join(raw,  by = c("study_num","substudy_num")) %>%
  filter(is.finite(beta), is.finite(vi))

# Species factor if present
if ("g_sp" %in% names(all)) all$g_sp <- as.factor(all$g_sp)

# Add vi_use for stable fitting/bias plots
all <- with_vi_use(all)

# Optional phylogeny structure (OFF by default)
R_phy <- NULL
if (USE_PHYLO) {
  if (!requireNamespace("ape", quietly=TRUE)) stop("Install 'ape' for phylogeny.")
  tree_path <- locate_file(c(TREE_FILE, file.path("output", TREE_FILE), file.path("data", TREE_FILE)))
  tree <- ape::read.tree(tree_path)
  keep_sp <- intersect(tree$tip.label, unique(all$g_sp))
  tree <- ape::keep.tip(tree, keep_sp)
  phylo_vcv <- ape::vcv(tree, corr=TRUE)
  keep_sp <- intersect(rownames(phylo_vcv), unique(all$g_sp))
  phylo_vcv <- phylo_vcv[keep_sp, keep_sp, drop=FALSE]
  all <- all %>% filter(g_sp %in% keep_sp)
  R_phy <- list(g_sp = phylo_vcv)
}

# =============================== Run ==========================================
# Baseline β (no covariates)
run_core_beta(dat = all,
              label = "baseline_beta",
              R_arg = NULL)      # switch to R_arg = R_phy if USE_PHYLO=TRUE

message("Core β diagnostics complete. See: ", ROOT_FIG, " / ", ROOT_TAB)