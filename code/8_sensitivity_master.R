# Core β Meta-Analysis — Publication-bias Sensitivity (asinh scale)
# R port with Leave-One-Out (LOO) metric
# Author: Stier Lab (ACS) | R refactor by ChatGPT
# Date: 2025-09-30

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(jsonlite)
  library(here)          # << add this
})


# Optional (preferred) packages
HAVE_METAFOR <- requireNamespace("metafor", quietly = TRUE)
HAVE_COWPLOT <- requireNamespace("cowplot", quietly = TRUE)
HAVE_XTABLE  <- requireNamespace("xtable",  quietly = TRUE)

# ---------------------------- Global plotting theme ---------------------------
theme_set(theme_minimal(base_size = 13))
update_geom_defaults("point", list(size = 2.5, alpha = 0.7, stroke = 0.5))
theme_core <- theme(
  axis.title  = element_text(size = 14),
  axis.text   = element_text(size = 12),
  plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.title= element_blank()
)

# ---------------------------- PLOTTING GLOBALS --------------------------------
set.seed(2025)
USE_JITTER    <- TRUE
JITTER_FRAC_X <- 0.010  # ~1% of x-range
JITTER_FRAC_Y <- 0.015  # ~1.5% of y-range

jitter_xy <- function(x, y, frac_x = JITTER_FRAC_X, frac_y = JITTER_FRAC_Y) {
  x <- as.numeric(x); y <- as.numeric(y)
  xrng <- diff(stats::quantile(x, c(.01, .99), na.rm = TRUE))
  yrng <- diff(stats::quantile(y, c(.01, .99), na.rm = TRUE))
  dx <- stats::rnorm(length(x), 0, pmax(1e-12, frac_x * xrng))
  dy <- stats::rnorm(length(y), 0, pmax(1e-12, frac_y * yrng))
  list(x = x + dx, y = y + dy)
}

# How to display the SE axis in funnels/diagnostics: 'log' | 'asinh' | 'linear'
SE_AXIS_SCALE <- "asinh"

nice_log_ticks <- function(ymin, ymax, max_ticks = 4) {
  ymin <- max(as.numeric(ymin), 1e-12)
  ymax <- max(as.numeric(ymax), ymin * (1 + 1e-12))
  lo <- floor(log10(ymin)); hi <- ceiling(log10(ymax))
  exps <- lo:hi
  if (length(exps) <= max_ticks) return(10^exps)
  step <- ceiling(length(exps) / max_ticks)
  kept <- unique(sort(c(lo, exps[seq(1, length(exps), by = step)], hi)))
  while (length(kept) > max_ticks) {
    mid <- 1 + floor((length(kept) - 2) / 2)
    kept <- kept[-mid]
  }
  10^kept
}

format_se_axis <- function(se_vals) {
  if (SE_AXIS_SCALE == "asinh") {
    ticks <- nice_log_ticks(min(se_vals, na.rm = TRUE), max(se_vals, na.rm = TRUE), max_ticks = 4)
    list(breaks = asinh(ticks), labels = sprintf("%g", ticks))
  } else if (SE_AXIS_SCALE == "log") {
    # We'll plot on original scale and set coord_trans, so labels are raw
    ticks <- nice_log_ticks(min(se_vals, na.rm = TRUE), max(se_vals, na.rm = TRUE), max_ticks = 5)
    list(breaks = ticks, labels = sprintf("%g", ticks))
  } else {
    NULL
  }
}


# ---- robust script directory detection (works in Rscript, RStudio, or interactive) ----
here_dir <- script_dir()

# ---- Prefer loading merged data by sourcing the pipeline step ----------------
load_from_phylogeny_loading <- function(here_dir) {
  # plausible locations/names for the data-prep script
  candidates <- c(
    file.path(here_dir, "1_data_phylogeny_loading.R"),
    file.path(here_dir, "R", "1_data_phylogeny_loading.R"),
    Sys.glob(file.path(here_dir, "*data*phylogeny*loading*.R"))
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (!length(candidates)) return(NULL)
  
  for (p in candidates) {
    env <- new.env(parent = emptyenv())
    ok <- tryCatch({ sys.source(p, envir = env); TRUE }, error = function(e) FALSE)
    if (!ok) next
    
    # prefer data frames; allow character path as a fallback
    for (nm in c("all_dat2", "all_dat", "alldat2", "alldat")) {
      if (exists(nm, envir = env, inherits = FALSE)) {
        obj <- get(nm, envir = env)
        if (is.data.frame(obj)) {
          df <- clean_names(tibble::as_tibble(obj))
          return(list(df = df, source = sprintf("<%s::%s>", basename(p), nm)))
        }
        if (is.character(obj) && length(obj) == 1 && file.exists(obj)) {
          df <- readr::read_csv(obj, show_col_types = FALSE) |> clean_names()
          return(list(df = df, source = normalizePath(obj, winslash = "/", mustWork = FALSE)))
        }
      }
    }
  }
  NULL
}




# ---------------------------- IO ----------------------------------------------
here_dir <- script_dir()

dir.create(file.path(here_dir, "figs",   "sensitivity_core_beta", "baseline_beta_pub"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(here_dir, "tables", "sensitivity_core_beta", "baseline_beta_pub"), recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(here_dir, "figs",   "sensitivity_core_beta", "baseline_beta_pub")
TAB_DIR <- file.path(here_dir, "tables", "sensitivity_core_beta", "baseline_beta_pub")

# Fallback globber (kept, but used only if sourcing fails)
latest_any <- function(patterns) {
  cand <- character(0)
  for (pat in patterns) {
    cand <- c(cand,
              Sys.glob(file.path(here_dir, pat)),
              Sys.glob(file.path(here_dir, "data", pat)),
              Sys.glob(file.path(here_dir, "output", pat)))
  }
  if (!length(cand)) return(NULL)
  cand[which.max(file.info(cand)$mtime)]
}

# ---------------------------- Load & Prepare ----------------------------------
# 1) Try sourcing the pipeline step to get all_dat2 in-memory
from_pipe <- load_from_phylogeny_loading(here_dir)

if (!is.null(from_pipe)) {
  alldat <- from_pipe$df
  ALLDAT_PATH <- from_pipe$source
  message(sprintf("[IO] Using merged data from pipeline: %s", ALLDAT_PATH))
} else {
  # 2) Fallback to CSVs in ., data/, or output/
  ALLDAT_PATH <- latest_any(c("all_dat2.csv","all_dat2*.csv","merged-covariates*.csv","merged_covariates*.csv"))
  if (is.null(ALLDAT_PATH)) {
    stop("Could not find merged data: either source 1_data_phylogeny_loading.R to create all_dat2, ",
         "or place a CSV matching all_dat2*.csv / merged-covariates*.csv in ., data/, or output/.",
         call. = FALSE)
  }
  message(sprintf("[IO] Using merged data file: %s", basename(ALLDAT_PATH)))
  alldat <- readr::read_csv(ALLDAT_PATH, show_col_types = FALSE) |> clean_names()
}

# ---------------------------- Utilities ---------------------------------------
clean_names <- function(df) {
  names(df) <- names(df) |>
    stringr::str_trim() |>
    str_replace_all(" ", "_") |>
    str_replace_all("-", "_") |>
    str_replace_all("\\.", "_") |>
    tolower()
  df
}

normalize_id <- function(v) {
  vapply(v, function(x) {
    if (is.na(x)) return(NA_character_)
    if (is.numeric(x)) {
      if (isTRUE(all.equal(x, round(x)))) return(as.character(as.integer(x)))
      return(as.character(x))
    }
    as.character(x)
  }, FUN.VALUE = character(1))
}

# ---------------------------- Load & Prepare ----------------------------------
alldat <- readr::read_csv(ALLDAT_PATH, show_col_types = FALSE) |> clean_names()

for (col in c("substudy_num", "study_num")) {
  if (col %in% names(alldat)) alldat[[col]] <- normalize_id(alldat[[col]])
}

required <- c("substudy_num","beta_hat","beta_variance_nls2")
if (!all(required %in% names(alldat))) {
  stop(sprintf("Missing required columns in %s: need %s, found %s",
               basename(ALLDAT_PATH), paste(required, collapse=", "), paste(names(alldat), collapse=", ")), call. = FALSE)
}

# Unit conversion and asinh transform + delta-method variance
alldat <- alldat |>
  mutate(
    beta_cm2 = beta_hat * 1e4,
    var_cm2  = beta_variance_nls2 * 1e8,
    beta     = asinh(beta_cm2),
    vi       = var_cm2 / (1 + beta_cm2^2)
  )

keep_cols <- c("substudy_num","beta","vi")
if ("study_num" %in% names(alldat)) keep_cols <- c(keep_cols, "study_num")
if ("g_sp"      %in% names(alldat)) keep_cols <- c(keep_cols, "g_sp")

dat <- alldat |>
  dplyr::select(dplyr::all_of(keep_cols)) |>
  distinct(substudy_num, .keep_all = TRUE) |>
  filter(is.finite(beta), is.finite(vi), vi > 0) |>
  mutate(sei = sqrt(vi)) |>
  arrange(substudy_num)

message(sprintf("[LOAD] %d usable substudies after QC", nrow(dat)))

# ---------------------------- Random-effects summary --------------------------
re_summary <- function(beta, vi) {
  beta <- as.numeric(beta); vi <- as.numeric(vi)
  ok <- is.finite(beta) & is.finite(vi) & (vi > 0)
  beta <- beta[ok]; vi <- vi[ok]
  if (HAVE_METAFOR) {
    # REML intercept-only meta-analysis
    m <- try(metafor::rma.uni(yi = beta, vi = vi, method = "REML"), silent = TRUE)
    if (!inherits(m, "try-error")) {
      mu <- as.numeric(m$b[1]); se <- m$se; tau2 <- m$tau2
      QE <- m$QE; QEp <- m$QEp; k <- m$k
      ci <- mu + c(-1,1) * 1.96 * se
      I2 <- max(0, (QE - (k - 1)) / QE) * 100
      return(list(mu = mu, se = se, ci_lb = ci[1], ci_ub = ci[2],
                  tau2 = tau2, Q = QE, Q_df = k - 1, Q_p = QEp,
                  I2 = I2, k = k))
    }
  }
  # DL fallback
  w <- 1/vi
  ybar <- sum(w * beta) / sum(w)
  Q <- sum(w * (beta - ybar)^2)
  df <- length(beta) - 1
  c  <- sum(w) - sum(w^2)/sum(w)
  tau2 <- ifelse(c > 0, pmax(0, (Q - df)/c), 0)
  w_star <- 1/(vi + tau2)
  mu <- sum(w_star * beta)/sum(w_star)
  se <- sqrt(1/sum(w_star))
  I2 <- if (Q > 0) pmax(0, (Q - df)/Q) * 100 else 0
  list(mu = mu, se = se, ci_lb = mu - 1.96*se, ci_ub = mu + 1.96*se,
       tau2 = tau2, Q = Q, Q_df = df, Q_p = NA_real_,
       I2 = I2, k = length(beta))
}

RE <- re_summary(dat$beta, dat$vi)

# ---------------------------- Egger test & plot -------------------------------
egger_fit <- function(df, zoom95 = FALSE) {
  d <- df |> filter(is.finite(sei), sei > 0)
  if (zoom95) {
    cap <- stats::quantile(d$sei, 0.95, na.rm = TRUE)
    d <- d |> filter(sei <= cap)
  }
  d <- d |> mutate(zi = beta / sei, prec = 1 / sei)
  intercept <- slope <- tval <- pval <- NA_real_
  if (nrow(d) >= 3) {
    m <- stats::lm(zi ~ prec, data = d)
    co <- summary(m)$coefficients
    intercept <- unname(co["(Intercept)", "Estimate"])
    slope     <- unname(co["prec", "Estimate"])
    tval      <- unname(co["(Intercept)", "t value"])
    pval      <- unname(co["(Intercept)", "Pr(>|t|)"])
  }
  list(intercept = intercept, slope = slope, t = tval, p = pval, n = nrow(d),
       scope = ifelse(zoom95, "zoom95", "full"))
}

egger_plot <- function(df, title, out_png, zoom95 = FALSE) {
  d <- df |> filter(is.finite(sei), sei > 0)
  if (zoom95) {
    cap <- stats::quantile(d$sei, 0.95, na.rm = TRUE)
    d <- d |> filter(sei <= cap)
  }
  d <- d |> mutate(zi = beta / sei, prec = 1 / sei)
  
  xcap <- stats::quantile(d$prec, 0.99, na.rm = TRUE)
  d_plot <- d |> filter(prec <= xcap)
  trimmed_n <- nrow(d) - nrow(d_plot)
  
  # model & prediction band
  intercept <- se_intercept <- pval <- NA_real_
  pred_df <- NULL
  if (nrow(d) >= 3) {
    m <- stats::lm(zi ~ prec, data = d)
    intercept <- coef(m)[1]; se_intercept <- sqrt(diag(vcov(m)))[1]
    pval <- summary(m)$coefficients[1,4]
    xg <- tibble(prec = seq(min(d_plot$prec), max(d_plot$prec), length.out = 200))
    pr <- predict(m, newdata = xg, interval = "prediction", level = 0.95)
    pred_df <- cbind(xg, as.data.frame(pr))
  }
  
  p <- ggplot() +
    { if (trimmed_n > 0)
      geom_point(data = d |> filter(prec > xcap), aes(prec, zi),
                 color = "grey70", stroke = 0.4) } +
    geom_point(data = d_plot, aes(prec, zi), shape = 21, fill = "#1f77b4", color = "black") +
    { if (!is.null(pred_df))
      list(
        geom_line(data = pred_df, aes(prec, fit), linetype = "dashed", linewidth = 1),
        geom_ribbon(data = pred_df, aes(prec, ymin = lwr, ymax = upr), alpha = 0.25)
      ) } +
    { if (is.finite(intercept) && is.finite(se_intercept))
      geom_errorbar(aes(x = 0, ymin = intercept - 1.96*se_intercept, ymax = intercept + 1.96*se_intercept),
                    width = 0.02) } +
    { if (is.finite(intercept))
      geom_point(aes(x = 0, y = intercept), shape = 22, fill = "white", color = "black", size = 2.5) } +
    labs(x = "Precision (1/SE)", y = "Standardized effect (asinh(β)/SE)",
         title = if (is.finite(pval))
           sprintf("%s\nintercept = %.2f (p = %.2e); n = %d",
                   if (zoom95) "Egger regression (zoom95)" else "Egger regression (full)",
                   intercept, pval, nrow(d))
         else (if (zoom95) "Egger regression (zoom95)" else "Egger regression (full)")
    ) +
    coord_cartesian(xlim = c(0, xcap * 1.02)) +
    theme_core +
    theme(panel.grid.minor = element_blank())
  
  if (trimmed_n > 0) {
    p <- p + annotate("text", x = xcap*1.02, y = min(d_plot$zi), hjust = 1, vjust = -0.2,
                      label = sprintf("trimmed %d at right 1%% of precision", trimmed_n), size = 3)
  }
  ggsave(out_png, p, width = 8, height = 6, dpi = 450)
  ggsave(sub("\\.png$", ".pdf", out_png), p, width = 8, height = 6, dpi = 450)
  list(intercept = intercept,
       slope = if (!is.na(pval)) unname(coef(m)[["prec"]]) else NA_real_,
       t = NA_real_, p = pval, n = nrow(d), scope = ifelse(zoom95, "zoom95","full"))
}

egger_full <- egger_plot(dat, "Egger regression (full; see table for exact p)",
                         file.path(FIG_DIR, "egger_full.png"), zoom95 = FALSE)
egger_95   <- egger_plot(dat, "Egger regression (zoom95; see table for exact p)",
                         file.path(FIG_DIR, "egger_zoom95.png"), zoom95 = TRUE)

# ---------------------------- Funnel helpers ----------------------------------
add_funnel_cones <- function(p, mu, se_vals, y_plot_vals) {
  # cones: mu ± 1.96*SE
  up  <- mu + 1.96 * se_vals
  low <- mu - 1.96 * se_vals
  dfc <- tibble(x = c(low, rev(up)), y = c(y_plot_vals, rev(y_plot_vals)))
  p + geom_path(data = tibble(x = up,  y = y_plot_vals), aes(x, y), linetype = "dotted") +
    geom_path(data = tibble(x = low, y = y_plot_vals), aes(x, y), linetype = "dotted")
}

# Contour shading (horizontal ribbons): build many thin rectangles
contour_ribbons_df <- function(mu, se_grid, z_levels = c(1.64, 1.96, 2.58), se_to_plot) {
  # se_to_plot is either identity(se) or asinh(se), precomputed
  purrr::map_dfr(z_levels, function(z) {
    tibble(
      z = z,
      ymin = head(se_to_plot, -1),
      ymax = tail(se_to_plot, -1),
      xmin = mu - z * head(se_grid, -1),
      xmax = mu + z * head(se_grid, -1)
    )
  })
}

# ---------------------------- Funnel Plot Function ----------------------------
save_funnel <- function(df, title, path, cap_quantile = NULL, add_hlines = TRUE, note = NULL) {
  d <- df
  if (!is.null(cap_quantile)) {
    se_cap <- stats::quantile(d$sei, cap_quantile, na.rm = TRUE)
    dropped <- d |> filter(sei > se_cap)
    d <- d |> filter(sei <= se_cap)
  } else {
    dropped <- d[0,]
  }
  d <- d |> mutate(sei_t = if (SE_AXIS_SCALE == "asinh") asinh(sei) else sei)
  
  # Points (with optional jitter)
  xy <- if (USE_JITTER) jitter_xy(d$beta, d$sei_t) else list(x = d$beta, y = d$sei_t)
  base <- tibble(beta = xy$x, sei_t = xy$y)
  
  dropped2 <- if (nrow(dropped) > 0) dplyr::mutate(dropped, sei_t = if (SE_AXIS_SCALE == "asinh") asinh(sei) else sei) else dropped
  p <- ggplot() +
    { if (nrow(dropped2) > 0) geom_point(data = dropped2, aes(beta, sei_t), color = "grey70", shape = 21, fill = "grey80") } +
    geom_point(data = base, aes(beta, sei_t), shape = 21, fill = "#1f77b4", color = "black") +
    annotate("rect", xmin = RE$ci_lb, xmax = RE$ci_ub,
             ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "orange") +
    geom_vline(xintercept = RE$mu, linetype = "dashed", color = "grey30", linewidth = 1) +
    labs(x = "asinh(β)", y = "Standard Error",
         title = sprintf("%s\n(n=%d, \u03bc=%.3f, 95%% CI [%.3f, %.3f])",
                         title, nrow(df), RE$mu, RE$ci_lb, RE$ci_ub)) +
    coord_cartesian(xlim = c(-8, 10)) +
    theme_core +
    theme(panel.grid.minor = element_blank())
  
  # SE percentiles
  if (add_hlines) {
    y95 <- stats::quantile(df$sei, 0.95, na.rm = TRUE)
    y99 <- stats::quantile(df$sei, 0.99, na.rm = TRUE)
    y95_t <- if (SE_AXIS_SCALE == "asinh") asinh(y95) else y95
    y99_t <- if (SE_AXIS_SCALE == "asinh") asinh(y99) else y99
    p <- p + geom_hline(yintercept = y95_t, linetype = "dashed", color = "grey50", linewidth = 1) +
      geom_hline(yintercept = y99_t, linetype = "dotted", color = "grey50", linewidth = 1)
  }
  
  # Cones
  se_vals  <- seq(min(d$sei, na.rm = TRUE), max(d$sei, na.rm = TRUE), length.out = 200)
  y_plot_v <- if (SE_AXIS_SCALE == "asinh") asinh(se_vals) else se_vals
  p <- add_funnel_cones(p, RE$mu, se_vals, y_plot_v)
  
  # Reverse y-axis (funnel style)
  p <- p + scale_y_continuous(trans = "reverse")
  
  # y ticks (raw SE labels even if asinh)
  se_axis <- format_se_axis(d$sei)
  if (!is.null(se_axis)) {
    p <- p + scale_y_continuous(trans = "reverse",
                                breaks = se_axis$breaks, labels = se_axis$labels)
  }
  
  if (!is.null(cap_quantile)) {
    p <- p + annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
                      label = sprintf("Trimmed %d SE>=%d%%", nrow(dropped), as.integer(cap_quantile*100)),
                      color = "firebrick", size = 3)
  }
  if (!is.null(note)) {
    p <- p + annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
                      label = note, size = 3)
  }
  
  ggsave(path, p, width = 7.5, height = 6, dpi = 450)
  ggsave(sub("\\.png$", ".pdf", path), p, width = 7.5, height = 6, dpi = 450)
  invisible(p)
}

save_contour_funnel <- function(df, title, path, z_levels = c(1.64, 1.96, 2.58), cap_quantile = NULL) {
  d <- df
  if (!is.null(cap_quantile)) {
    se_cap <- stats::quantile(d$sei, cap_quantile, na.rm = TRUE)
    d <- d |> filter(sei <= se_cap)
  }
  d <- d |> mutate(sei_t = if (SE_AXIS_SCALE == "asinh") asinh(sei) else sei)
  
  # grid for ribbons
  q_cap <- if (is.null(cap_quantile)) 0.95 else min(0.95, cap_quantile)
  se_grid <- seq(min(d$sei, na.rm = TRUE), stats::quantile(d$sei, q_cap, na.rm = TRUE), length.out = 400)
  y_plot <- if (SE_AXIS_SCALE == "asinh") asinh(se_grid) else se_grid
  rib <- contour_ribbons_df(RE$mu, se_grid, z_levels = sort(z_levels), se_to_plot = y_plot) |>
    mutate(lbl = dplyr::case_when(
      abs(z - 1.64) < 1e-6 ~ "p<0.10",
      abs(z - 1.96) < 1e-6 ~ "p<0.05",
      abs(z - 2.58) < 1e-6 ~ "p<0.01",
      TRUE ~ sprintf("z=%.2f", z)
    ))
  
  # jittered points
  dxy <- if (USE_JITTER) jitter_xy(d$beta, d$sei_t) else list(x = d$beta, y = d$sei_t)
  base_pts <- tibble(beta = dxy$x, sei_t = dxy$y)
  
  p <- ggplot() +
    geom_rect(data = rib, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = lbl),
              alpha = 0.25, color = NA) +
    scale_fill_manual(values = c("p<0.10" = "grey80", "p<0.05" = "lightsteelblue3", "p<0.01" = "lightblue"),
                      drop = FALSE) +
    geom_point(data = base_pts, aes(beta, sei_t), shape = 21, fill = "#1f77b4", color = "black") +
    annotate("rect", xmin = RE$ci_lb, xmax = RE$ci_ub,
             ymin = -Inf, ymax = Inf, alpha = 0.20, fill = "orange") +
    geom_vline(xintercept = RE$mu, linetype = "dashed") +
    labs(x = "asinh(β)", y = "Standard Error",
         title = sprintf("%s\n(n=%d, \u03bc=%.3f, 95%% CI [%.3f, %.3f])",
                         title, nrow(df), RE$mu, RE$ci_lb, RE$ci_ub)) +
    coord_cartesian(xlim = c(-8, 10)) +
    theme_core +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  
  # reverse y and raw SE ticks
  p <- p + scale_y_continuous(trans = "reverse")
  se_axis <- format_se_axis(d$sei)
  if (!is.null(se_axis)) {
    p <- p + scale_y_continuous(trans = "reverse",
                                breaks = se_axis$breaks, labels = se_axis$labels)
  }
  
  ggsave(path, p, width = 7.5, height = 6, dpi = 450)
  ggsave(sub("\\.png$", ".pdf", path), p, width = 7.5, height = 6, dpi = 450)
  invisible(p)
}

# ---------------------------- Make plots --------------------------------------
save_funnel(dat, "A) Funnel (All points, asinh β)",
            file.path(FIG_DIR, "funnel_A_all.png"),
            cap_quantile = NULL, add_hlines = TRUE,
            note = "Most SE values are small; see zoomed panels.")

save_contour_funnel(dat, "A0) Contour-enhanced funnel (asinh β)",
                    file.path(FIG_DIR, "funnel_A0_contour.png"), cap_quantile = 0.95)

# minus the single largest-SE point
idx_max <- which.max(dat$sei)
dat_minus1 <- dat[-idx_max, , drop = FALSE]
save_funnel(dat_minus1, "A2) Funnel (All points minus largest-SE outlier)",
            file.path(FIG_DIR, "funnel_A2_minus_top.png"),
            cap_quantile = NULL, add_hlines = TRUE,
            note = "Only the single largest-SE point removed.")

save_funnel(dat, "B) Funnel (Zoom95, asinh β)",
            file.path(FIG_DIR, "funnel_B_zoom95.png"), cap_quantile = 0.95, add_hlines = TRUE)

save_funnel(dat, "C) Funnel (Capped 99.5% SE, asinh β)",
            file.path(FIG_DIR, "funnel_C_cap995.png"), cap_quantile = 0.995, add_hlines = TRUE)

# Scatter reference
scatter_ref <- function(df, out) {
  sei_t <- if (SE_AXIS_SCALE == "asinh") asinh(df$sei) else df$sei
  xy <- if (USE_JITTER) jitter_xy(df$beta, sei_t) else list(x = df$beta, y = sei_t)
  y95 <- stats::quantile(df$sei, 0.95, na.rm = TRUE)
  y99 <- stats::quantile(df$sei, 0.99, na.rm = TRUE)
  y95_t <- if (SE_AXIS_SCALE == "asinh") asinh(y95) else y95
  y99_t <- if (SE_AXIS_SCALE == "asinh") asinh(y99) else y99
  
  p <- ggplot(tibble(beta = xy$x, sei_t = xy$y), aes(beta, sei_t)) +
    geom_point(shape = 21, fill = "#1f77b4", color = "black") +
    geom_hline(yintercept = y95_t, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = y99_t, linetype = "dotted", color = "grey50") +
    scale_y_continuous(trans = "reverse") +
    { se_axis <- format_se_axis(df$sei);
    if (!is.null(se_axis)) scale_y_continuous(trans = "reverse", breaks = se_axis$breaks, labels = se_axis$labels) } +
    labs(x = "asinh(β)", y = "Standard Error", title = "D) asinh(β) vs SE (95/99 cutoffs)") +
    theme_core + coord_cartesian(xlim = c(-8, 10))
  ggsave(out, p, width = 7.5, height = 6, dpi = 450)
  ggsave(sub("\\.png$", ".pdf", out), p, width = 7.5, height = 6, dpi = 450)
}
scatter_ref(dat, file.path(FIG_DIR, "funnel_D_scatter.png"))

# ---------------------------- Begg rank-correlation test ----------------------
begg_stats <- list(n = nrow(dat), tau = NA_real_, p = NA_real_)
if (nrow(dat) >= 10) {
  # Kendall correlation between |beta| and SE
  kt <- suppressWarnings(cor.test(abs(dat$beta), dat$sei, method = "kendall", exact = FALSE))
  begg_stats$tau <- unname(kt$estimate)
  begg_stats$p   <- unname(kt$p.value)
}

# ---------------------------- Influence / Sensitivity -------------------------
re_dl <- function(beta, vi) {
  beta <- as.numeric(beta); vi <- as.numeric(vi)
  ok <- is.finite(beta) & is.finite(vi) & (vi > 0)
  beta <- beta[ok]; vi <- vi[ok]
  w <- 1/vi
  ybar <- sum(w * beta) / sum(w)
  Q <- sum(w * (beta - ybar)^2)
  df <- length(beta) - 1
  c  <- sum(w) - sum(w^2)/sum(w)
  tau2 <- ifelse(c > 0, pmax(0, (Q - df)/c), 0)
  sum(beta/(vi + tau2)) / sum(1/(vi + tau2))
}

mu_full <- RE$mu

# LOO by substudy (remove ALL rows with that substudy_num—already 1/row)
loo_df <- dat |>
  group_by(substudy_num) |>
  group_modify(function(g, .) {
    mask <- dat$substudy_num != g$substudy_num[1]
    mu_i <- re_dl(dat$beta[mask], dat$vi[mask])
    tibble(loo_intercept_asinh = mu_i)
  }) |>
  ungroup()

if ("study_num" %in% names(dat)) loo_df$study_num <- dat$study_num
loo_df$delta <- loo_df$loo_intercept_asinh - mu_full

# Plot LOO series
p_loo <- ggplot(loo_df, aes(seq_along(loo_intercept_asinh), loo_intercept_asinh)) +
  geom_point(color = "#1f77b4") +
  geom_hline(yintercept = mu_full, linetype = "dashed", color = "grey30") +
  labs(x = "Left-out substudy index", y = "Intercept (asinh β)",
       title = "Leave-one-out meta-analytic intercept by substudy (DerSimonian–Laird)") +
  theme_core
ggsave(file.path(FIG_DIR, "leave_one_out.png"), p_loo, width = 12, height = 4.8, dpi = 450)
ggsave(file.path(FIG_DIR, "leave_one_out.pdf"), p_loo, width = 12, height = 4.8, dpi = 450)

# LOO by paper (study_num), if present
if ("study_num" %in% names(dat)) {
  loo_study <- dat |>
    group_by(study_num) |>
    group_map(function(g, key) {
      mask <- dat$study_num != key$study_num
      mu_i <- re_dl(dat$beta[mask], dat$vi[mask])
      tibble(study_num = key$study_num, k_removed = nrow(g),
             loo_intercept_asinh = mu_i, delta = mu_i - mu_full)
    }) |>
    list_rbind() |>
    arrange(desc(abs(delta)))
}

# LOO by species (g_sp), if present
if ("g_sp" %in% names(dat)) {
  loo_sp <- dat |>
    group_by(g_sp) |>
    group_map(function(g, key) {
      mask <- dat$g_sp != key$g_sp
      mu_i <- re_dl(dat$beta[mask], dat$vi[mask])
      tibble(g_sp = key$g_sp, k_removed = nrow(g),
             loo_intercept_asinh = mu_i, delta = mu_i - mu_full)
    }) |>
    list_rbind() |>
    arrange(desc(abs(delta)))
}

# Worst-case "file drawer": add m nulls at low precision
null_add_sensitivity <- function(df, m_values = c(1,3,5,10), se_quantile = 0.95) {
  se0 <- as.numeric(stats::quantile(df$sei, se_quantile, na.rm = TRUE))
  vi0 <- se0^2
  purrr::map_dfr(m_values, function(m) {
    mu_i <- re_dl(c(df$beta, rep(0, m)), c(df$vi, rep(vi0, m)))
    tibble(m_added_nulls = m, assumed_se = se0, loo_intercept_asinh = mu_i,
           delta = mu_i - mu_full)
  })
}
null_grid <- null_add_sensitivity(dat)

# ---------------------------- Unified summary table ---------------------------
build_summary_df <- function() {
  top15 <- loo_df |>
    mutate(idx = row_number()) |>
    slice_max(order_by = abs(delta), n = 15, with_ties = FALSE)
  
  summary_row <- tibble(
    label = "baseline_beta_asinh",
    k = RE$k,
    intercept_asinh = RE$mu,
    intercept_se_asinh = RE$se,
    ci_lb_asinh = RE$ci_lb,
    ci_ub_asinh = RE$ci_ub,
    mean_beta_raw = sinh(RE$mu),
    ci_lb_beta_raw = sinh(RE$ci_lb),
    ci_ub_beta_raw = sinh(RE$ci_ub),
    QE = RE$Q,
    QEp = RE$Q_p,
    I2_overall = RE$I2,
    tau2_components = RE$tau2,
    tau2_names = "overall",
    egger_full_intercept = egger_full$intercept,
    egger_full_slope     = egger_full$slope,
    egger_full_t         = egger_full$t,
    egger_full_p         = egger_full$p,
    egger_full_n         = egger_full$n,
    egger_zoom95_intercept = egger_95$intercept,
    egger_zoom95_slope     = egger_95$slope,
    egger_zoom95_t         = egger_95$t,
    egger_zoom95_p         = egger_95$p,
    egger_zoom95_n         = egger_95$n,
    loo_mu_min = min(loo_df$loo_intercept_asinh),
    loo_mu_max = max(loo_df$loo_intercept_asinh),
    loo_mu_sd  = stats::sd(loo_df$loo_intercept_asinh),
    loo_delta_max_abs = max(abs(loo_df$delta)),
    begg_tau = begg_stats$tau,
    begg_p   = begg_stats$p,
    delta_with_5_nulls = null_grid |> filter(m_added_nulls == 5) |> pull(delta)
  )
  
  bind_rows(
    summary_row,
    mutate(loo_df, label = "leave_one_out"),
    mutate(top15,  label = "leave_one_out_top15") |> select(-idx),
    mutate(null_grid, label = "null_add_sensitivity"),
    tibble(label = "begg_rank_correlation", n = nrow(dat), tau = begg_stats$tau, p = begg_stats$p)
  )
}

summary_df <- build_summary_df()

summary_path <- file.path(TAB_DIR, "AppendixD_summary.csv")
readr::write_csv(summary_df, summary_path, na = "")

# Optional LaTeX using xtable (fallback: skip)
if (HAVE_XTABLE) {
  texfile <- file.path(TAB_DIR, "AppendixD_summary.tex")
  sink(texfile); on.exit(sink(), add = TRUE)
  print(xtable::xtable(summary_df, digits = 4), include.rownames = FALSE)
  on.exit(NULL, add = FALSE)
}

# Detailed CSVs
readr::write_csv(loo_df, file.path(TAB_DIR, "leave_one_out.csv"))
top15 <- loo_df |> slice_max(order_by = abs(delta), n = 15, with_ties = FALSE)
readr::write_csv(top15, file.path(TAB_DIR, "leave_one_out_top15.csv"))

if (exists("loo_study")) readr::write_csv(loo_study, file.path(TAB_DIR, "leave_one_study.csv"))
if (exists("loo_sp"))    readr::write_csv(loo_sp,    file.path(TAB_DIR, "leave_one_species.csv"))
readr::write_csv(null_grid, file.path(TAB_DIR, "sensitivity_add_nulls.csv"))
readr::write_csv(tibble(n = nrow(dat), tau = begg_stats$tau, p = begg_stats$p),
                 file.path(TAB_DIR, "begg_rank_correlation.csv"))

# RUN LOG
log_path <- file.path(TAB_DIR, "RUN_LOG.txt")
log_list <- list(RE = RE, Egger_full = egger_full, Egger_zoom95 = egger_95,
                 Begg = begg_stats, Null_sensitivity = null_grid)
writeLines(c(
  "Core β Meta-Analysis — R run log",
  sprintf("Inputs: %s", ALLDAT_PATH),
  jsonlite::toJSON(log_list, auto_unbox = TRUE, pretty = TRUE)
), con = log_path)

message("Done.\nFigures: ", FIG_DIR, "\nSummary table: ", summary_path)