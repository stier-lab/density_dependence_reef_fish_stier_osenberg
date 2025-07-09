# File: 3_phylo_species_analysis.R
# Purpose: Estimate species-level beta effects, map onto phylogeny, assess phylogenetic signal,
#          compare phylogenetic vs. non-phylogenetic models, partition variance,
#          wrasse example, and plot species & family-level patterns.
# Author: Adrian Stier
# Date: 2025-07-09

# ----------------------------------------------------------------------------
# 1. Dependencies & Data Load
# ----------------------------------------------------------------------------
source(here::here("code", "0_libraries.R"))               # load all required packages
source(here::here("code", "1_data_phylogeny_loading.R"))  # loads all_dat2, pruned_tree, etc.
source(here::here("code", "2_beta_estimate.R"))           # computes betanls2_raw_cm, betanls2_asinh, etc.

# Ensure output directories exist
dir.create(here::here("figures"), showWarnings = FALSE)
dir.create(here::here("results"), showWarnings = FALSE)

# Define color palette
pal <- wesanderson::wes_palette("Zissou1", n = 100, type = "continuous")
# ----------------------------------------------------------------------------
# 2. Prepare species-level dataset
# ----------------------------------------------------------------------------
species_df <- all_dat2 %>%
  dplyr::select(
    g_sp,
    betanls2_raw_cm,
    betanls2_asinh,
    betanlsvar_asinh,
    study_num,
    substudy_num,
    family,
    predators,
    paired_pred,
    expt_obs_pairs,
    duration,
    mean_density,
    max_length_cm,
    max_weight_g
  ) %>%
  # Standardize a few species names
  dplyr::mutate(
    g_sp = stringr::str_replace_all(g_sp, c(
      "Sebastes_spp_NA"           = "Sebastes_mystinus",
      "Fundulus_heteroclitus"     = "Fundulus_heteroclitus_heteroclitus",
      "Diplodus_sargus"           = "Diplodus_sargus_sargus",
      "Elactinus_sp\\."           = "Elacatinus_evelynae"
    ))
  ) %>%
  dplyr::filter(!is.na(betanls2_asinh), !is.na(betanlsvar_asinh))

# ----------------------------------------------------------------------------
# 3. Species-level meta-analysis (no phylogeny)
# ----------------------------------------------------------------------------
sp_model <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ g_sp - 1,                   # species-specific estimates
  random = list(~1 | study_num/substudy_num),
  data   = species_df,
  method = "REML"
)

sp_coefs <- coef(summary(sp_model)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  dplyr::rename(
    est_asinh   = estimate,
    ci_lb_asinh = ci.lb,
    ci_ub_asinh = ci.ub
  ) %>%
  dplyr::mutate(
    g_sp     = stringr::str_remove(term, "^g_sp"),
    beta_est = sinh(est_asinh),
    ci.lb    = sinh(ci_lb_asinh),
    ci.ub    = sinh(ci_ub_asinh)
  ) %>%
  dplyr::select(g_sp, beta_est, ci.lb, ci.ub)

readr::write_csv(sp_coefs, here::here("results", "species_effects_no_phylo.csv"))

# ----------------------------------------------------------------------------
# 4. Map species estimates onto phylogeny
# ----------------------------------------------------------------------------
# Create data frame for plotting
beta_plot_df <- sp_coefs %>%
  dplyr::rename(mean = beta_est)

# Simple circular tree with colored points
p_phylo_simple <- ggtree::ggtree(pruned_tree, layout = "circular") %<+% beta_plot_df +
  ggtree::geom_tippoint(aes(color = mean), size = 3) +
  ggplot2::scale_color_viridis_c(trans = "asinh", name = "Beta") +
  ggtree::geom_tiplab(size = 2, offset = 1.5) +
  ggplot2::theme(legend.position = "right")

# Detailed circular tree: black outline + white-to-red fill
p_phylo <- ggtree::ggtree(pruned_tree, layout = "circular") %<+% beta_plot_df +
  ggtree::geom_tippoint(
    aes(fill = mean),
    shape = 21, size = 5, color = "black", stroke = 1
  ) +
  ggplot2::scale_fill_gradientn(
    colours = pal,
    name    = "Beta",
    trans   = "asinh",
    limits  = c(min(beta_plot_df$mean), max(beta_plot_df$mean)),
    breaks  = c(-10, 0, 10, 100, 1000, 10000),
    labels  = c("-10", "0", "10", "100", "1,000", "10,000"),
    guide   = guide_colorbar(
      frame.colour    = "black",
      frame.linewidth = 0.5,
      ticks.colour    = "black",
      ticks.linewidth = 0.5
    )
  ) +
  ggtree::geom_tiplab(size = 2, offset = 2) +
  ggtree::theme_tree2() +
  ggplot2::theme(
    axis.line   = element_blank(),
    panel.grid  = element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  )

# Save the detailed plot
ggplot2::ggsave(
  filename = here::here("figures", "phylo_beta_circular.png"),
  plot     = p_phylo,
  width    = 8, height = 8, dpi = 300
)

# ----------------------------------------------------------------------------
# 5. Phylogenetic signal tests
# ----------------------------------------------------------------------------

betavec <- beta_plot_df$mean
names(betavec) <- beta_plot_df$g_sp

# only keep the species that appear in both the tree and your estimates:
common_sp <- intersect(pruned_tree$tip.label, names(betavec))

# subset and reorder
betavec <- betavec[common_sp]
pruned_tree <- ape::drop.tip(pruned_tree,
                             setdiff(pruned_tree$tip.label, common_sp))

lambda_res <- phytools::phylosig(pruned_tree, betavec,
                                 method = "lambda", test = TRUE)
k_res      <- phytools::phylosig(pruned_tree, betavec,
                                 method = "K",      test = TRUE)


# ----------------------------------------------------------------------------
# 6. Compare models: no-phylo vs. phylo-random
# ----------------------------------------------------------------------------
# Build phylogenetic VCV matrix
phylo_vcv <- ape::vcv(pruned_tree, corr = TRUE)

# Only keep species that are in the VCV matrix
species_phylo <- species_df %>%
  dplyr::filter(g_sp %in% rownames(phylo_vcv)) %>%
  droplevels()

# (a) no‐phylo model on the same data subset
m_nosp <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~1,
  random = list(~1 | study_num/substudy_num),
  data   = species_phylo,
  method = "REML"
)

# (b) phylo‐random model
m_gsp <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~1,
  random = list(~1 | study_num/substudy_num, ~1 | g_sp),
  R      = list(g_sp = phylo_vcv),
  data   = species_phylo,
  method = "REML"
)
comp_aic <- AIC(m_nosp, m_gsp)
readr::write_csv(as.data.frame(comp_aic), here::here("results", "model_comparison_aic.csv"))

# Z-test for species variance
species_var <- m_gsp$sigma2[3]
species_se  <- sqrt(species_var)
z_species   <- species_var / species_se
p_species   <- 2 * pnorm(-abs(z_species))

tibble::tibble(
  Z_score = z_species,
  P_value = p_species
) %>%
  readr::write_csv(here::here("results", "species_variance_test.csv"))

# ----------------------------------------------------------------------------
# 7. Variance partitioning
# ----------------------------------------------------------------------------
sigma      <- m_gsp$sigma2
names(sigma) <- c("study", "substudy", "species")
total_var  <- sigma["substudy"] + sigma["species"]

partition <- tibble::tibble(
  level      = names(sigma),
  variance   = sigma,
  proportion = case_when(
    level == "substudy" ~ sigma["substudy"] / total_var,
    level == "species"  ~ sigma["species"]  / total_var,
    TRUE                ~ NA_real_
  )
)
readr::write_csv(partition, here::here("results", "variance_partitioning.csv"))

# ----------------------------------------------------------------------------
# 8. Wrasse example: within vs. among variance
# ----------------------------------------------------------------------------
wrasse_data <- species_df %>%
  dplyr::filter(g_sp %in% c("Thalassoma_bifasciatum", "Halichoeres_garnoti")) %>%
  dplyr::select(g_sp, study_num, substudy_num, betanls2_raw_cm)

species_means <- wrasse_data %>%
  dplyr::group_by(g_sp) %>%
  dplyr::summarise(species_mean = mean(betanls2_raw_cm, na.rm = TRUE))

wrasse_data <- wrasse_data %>%
  dplyr::left_join(species_means, by = "g_sp") %>%
  dplyr::mutate(deviation = betanls2_raw_cm - species_mean)

within_var <- var(wrasse_data$deviation, na.rm = TRUE)
among_var  <- var(species_means$species_mean, na.rm = TRUE)
total_w    <- var(wrasse_data$betanls2_raw_cm, na.rm = TRUE)

wrasse_part <- tibble::tibble(
  component = c("within_species", "among_species", "total"),
  variance  = c(within_var, among_var, total_w),
  percent   = c(within_var, among_var, total_w) / total_w * 100
)
readr::write_csv(wrasse_part, here::here("results", "wrasse_variance.csv"))

# ----------------------------------------------------------------------------
# 9. Plot species by family ordered β
# ----------------------------------------------------------------------------
all2 <- species_df %>%
  dplyr::filter(predators == "present") %>%
  droplevels() %>%
  dplyr::rename(beta_raw = betanls2_raw_cm)

sp_coefs_plot <- sp_coefs %>%
  dplyr::left_join(all_dat2 %>% dplyr::distinct(g_sp, family), by = "g_sp")

family_order  <- sp_coefs_plot %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(mean_beta = mean(beta_est)) %>%
  dplyr::arrange(desc(mean_beta)) %>%
  dplyr::pull(family)

species_order <- sp_coefs_plot %>%
  dplyr::mutate(family = factor(family, levels = family_order)) %>%
  dplyr::arrange(family, desc(beta_est)) %>%
  dplyr::pull(g_sp)

all2$g_sp         <- factor(all2$g_sp,         levels = species_order)
sp_coefs_plot$g_sp <- factor(sp_coefs_plot$g_sp, levels = species_order)

palette9 <- wesanderson::wes_palette("Darjeeling1", type = "continuous", n = 9)

p_sp <- ggplot2::ggplot() +
  ggplot2::geom_jitter(
    data    = all2,
    ggplot2::aes(x = beta_raw, y = g_sp, colour = family),
    width   = 0.2, height = 0, size = 2, alpha = 0.5
  ) +
  ggplot2::geom_errorbarh(
    data    = sp_coefs_plot,
    ggplot2::aes(y = g_sp, xmin = ci.lb, xmax = ci.ub, colour = family),
    height  = 0, size = 0.8
  ) +
  ggplot2::geom_point(
    data    = sp_coefs_plot,
    ggplot2::aes(x = beta_est, y = g_sp, colour = family),
    shape   = 17, size = 4
  ) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::scale_x_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
    labels = scales::comma
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) parse(text = paste0("italic(", gsub("_", "~", x), ")"))
  ) +
  ggplot2::scale_colour_manual(
    values = palette9,
    breaks = family_order,
    name   = "Family"
  ) +
  ggplot2::labs(x = "Density-Dependence (β)", y = NULL) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.text.y     = element_text(size = 10),
    legend.position = "right"
  )

ggplot2::ggsave(
  filename = here::here("figures", "species_by_family_beta.png"),
  plot     = p_sp,
  width    = 10, height = 8, dpi = 300
)

# End of 3_phylo_species_analysis.R