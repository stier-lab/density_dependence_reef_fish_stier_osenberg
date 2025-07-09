# =============================================================================
# File: 0_libraries.R
# Purpose: Centralized, minimal set of packages for the entire phylogenetic meta-analysis
# Author: Adrian Stier
# Date: 2025-07-09
# =============================================================================



# ----------------------------
# 1. Data wrangling & I/O
# ----------------------------
required_pkgs_tidy <- c("dplyr", "tidyr", "purrr", "broom","readr", "tibble")
for (pkg in required_pkgs_tidy) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 2. Project structure & strings
# ----------------------------
for (pkg in c("here", "stringr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 3. Phylogenetics
# ----------------------------
for (pkg in c("ape", "phytools", "ggtree")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 4. Meta-analysis & modeling
# ----------------------------
if (!requireNamespace("metafor", quietly = TRUE)) {
  stop("Package 'metafor' is required but not installed.")
}
library(metafor)

# ----------------------------
# 5. Visualization
# ----------------------------
required_pkgs_viz <- c("ggplot2", "GGally", "scales", "viridis", "wesanderson","fishtree")
for (pkg in required_pkgs_viz) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 6. Tables & reporting
# ----------------------------
if (!requireNamespace("gt", quietly = TRUE)) {
  stop("Package 'gt' is required but not installed.")
}
library(gt)

# ----------------------------
# 7. Reproducibility
# ----------------------------
# consistent random seed
set.seed(20250708)

# print session info for reproducibility
message("--- Session information for reproducibility ---")
print(sessionInfo())