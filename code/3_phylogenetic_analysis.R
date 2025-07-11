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


# bootstrap 95% CI for Blomberg's K

nsim   <- 500
boot_K <- numeric(nsim)

for(i in seq_len(nsim)) {
  sim_trait   <- rTraitCont(pruned_tree)
  boot_K[i]   <- phylosig(pruned_tree, sim_trait, method="K", test=FALSE)
}

ci_K <- quantile(boot_K, c(0.025, 0.975))
ci_K


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
library(dplyr)

wrasse_data <- all_dat2 %>%
  filter(g_sp %in% c("Thalassoma_bifasciatum", "Halichoeres_garnoti")) %>%
  select(g_sp, study_num, substudy_num, betanls2_raw_cm)


# 1) Compute species means and overall mean
species_means <- wrasse_data %>%
  group_by(g_sp) %>%
  summarise(mu_sp = mean(betanls2_raw_cm), n_sp = n(), .groups="drop")

grand_mu <- with(species_means, sum(mu_sp * n_sp) / sum(n_sp))

# 2) Among-species variance
among_var <- sum( (species_means$mu_sp - grand_mu)^2 * species_means$n_sp ) /
  (sum(species_means$n_sp) - 1)

# 3) Within-species variance (weighted average)
within_var <- wrasse_data %>%
  left_join(species_means, by="g_sp") %>%
  mutate(dev = betanls2_raw_cm - mu_sp) %>%
  summarise(w = sum(dev^2) / (n() - n_distinct(g_sp))) %>%
  pull(w)

# 4) Percent partition
partition <- tibble(
  component = c("within_species", "among_species"),
  variance  = c(within_var, among_var)
) %>%
  mutate(percent = variance / sum(variance) * 100)

partition

# ----------------------------------------------------------------------------
# 9. Plot species by family ordered β
# ----------------------------------------------------------------------------


# ──────────────────────────────────────────────────────────────────────
# 1. Prepare the data exactly as for the model (incl. study IDs!)
# ──────────────────────────────────────────────────────────────────────
all2 <- all_dat2 %>%
  filter(!is.na(betanls2_asinh),
         !is.na(betanlsvar_asinh),
         predators == "present") %>%     # correct level
  rename(beta_raw = betanls2_raw_cm) %>%
  select(g_sp, study_num, substudy_num, beta_raw, betanls2_asinh, betanlsvar_asinh) %>%
  droplevels()

# ──────────────────────────────────────────────────────────────────────
# 2. Fit the random‐effects species model (no intercept)
# ──────────────────────────────────────────────────────────────────────
m_sp_ni <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ g_sp - 1,
  random = list(~1 | study_num/substudy_num),
  data   = all2,
  method = "REML",
  test   = "t"
)

# ──────────────────────────────────────────────────────────────────────
# 3. Pull out which species actually went into that model
# ──────────────────────────────────────────────────────────────────────
model_species <- sub("^g_sp","", rownames(coef(summary(m_sp_ni))))

# ──────────────────────────────────────────────────────────────────────
# 4. Extract & back‐transform estimates + CIs
# ──────────────────────────────────────────────────────────────────────
sp_coefs <- coef(summary(m_sp_ni)) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  rename(
    est_asinh   = estimate,
    ci.lb_asinh = ci.lb,
    ci.ub_asinh = ci.ub
  ) %>%
  mutate(
    g_sp     = sub("^g_sp","",term),
    beta_est = sinh(est_asinh),
    ci.lb    = sinh(ci.lb_asinh),
    ci.ub    = sinh(ci.ub_asinh)
  ) %>%
  select(g_sp, beta_est, ci.lb, ci.ub)

# ──────────────────────────────────────────────────────────────────────
# 5. Restrict both to the modeled species
# ──────────────────────────────────────────────────────────────────────
all2     <- all2    %>% filter(g_sp %in% model_species) %>% droplevels()
sp_coefs <- sp_coefs %>% filter(g_sp %in% model_species)

# ──────────────────────────────────────────────────────────────────────
# 6. Add fish‐family for coloring
# ──────────────────────────────────────────────────────────────────────
family_df <- all_dat2 %>% distinct(g_sp, family)
all2     <- left_join(all2,     family_df, by = "g_sp")
sp_coefs <- left_join(sp_coefs, family_df, by = "g_sp")

# ──────────────────────────────────────────────────────────────────────
# 7. Order families by mean β, then species by β within each family
# ──────────────────────────────────────────────────────────────────────
family_order <- sp_coefs %>%
  group_by(family) %>%
  summarise(mean_beta = mean(beta_est)) %>%
  arrange(desc(mean_beta)) %>%
  pull(family)

species_order <- sp_coefs %>%
  mutate(family = factor(family, levels = family_order)) %>%
  arrange(family, desc(beta_est)) %>%
  pull(g_sp)

all2$g_sp     <- factor(all2$g_sp,     levels = species_order)
sp_coefs$g_sp <- factor(sp_coefs$g_sp, levels = species_order)
# ──────────────────────────────────────────────────────────────────────
# After step 7 (you already have `family_order` and `species_order`)
# ──────────────────────────────────────────────────────────────────────

# 7b) make sure `family` is a factor in the same order
all2$family    <- factor(all2$family,    levels = family_order)
sp_coefs$family <- factor(sp_coefs$family, levels = family_order)

# 8a) drop any non‐finite β’s so geom_point/geom_jitter won’t silently remove them
all2     <- all2    %>% filter(is.finite(beta_raw))
sp_coefs <- sp_coefs %>% filter(is.finite(beta_est), is.finite(ci.lb), is.finite(ci.ub))

# ──────────────────────────────────────────────────────────────────────
# After you’ve built all2, sp_coefs, family_order & species_order…
# ──────────────────────────────────────────────────────────────────────

library(ggplot2)
library(wesanderson)

# 1) Build Darjeeling1 ×9
palette_darjeeling9 <- wes_palette("Darjeeling1",
                                   type = "continuous",
                                   n    = 9)

# 2) Re-draw your plot with that palette
p_sp_darjeeling <- ggplot() +
  # raw sub‐study dots
  geom_jitter(
    data    = all2,
    aes(x = beta_raw, y = g_sp, colour = family),
    width   = 0.2, height = 0, size = 2, alpha = 0.5
  ) +
  # species‐level CIs
  geom_errorbarh(
    data    = sp_coefs,
    aes(y = g_sp, xmin = ci.lb, xmax = ci.ub, colour = family),
    height = 0, size = 0.8
  ) +
  # species‐level point estimates
  geom_point(
    data    = sp_coefs,
    aes(x = beta_est, y = g_sp, colour = family),
    shape = 17, size = 4
  ) +
  # vertical zero line
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  # asinh x‐axis
  scale_x_continuous(
    trans  = "asinh",
    breaks = c(-1000,-100,-10,-1,0,1,10,100,1000,10000),
    labels = scales::comma
  ) +
  # italic y‐axis labels (no underscores)
  scale_y_discrete(
    labels = function(x) {
      txt <- gsub("_","~", x)
      parse(text = paste0("italic(", txt, ")"))
    }
  ) +
  # Darjeeling1×9 manual palette, legend in same top→bottom order
  scale_colour_manual(
    values = palette_darjeeling9,
    breaks = family_order,
    name   = "Fish Family",
    guide  = guide_legend(reverse = TRUE)
  ) +
  labs(
    x = expression(
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
    y = "",
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y        = element_text(size = 10),
    axis.title.x       = element_text(size = 12, face = "bold"),
    legend.position    = "right",
    panel.grid.major.x = element_line(color = "gray90", linetype = "dotted")
  )

# 3) Render & save
print(p_sp_darjeeling)
ggsave("figures/Figure3_species_by_family_order_beta_darjeeling9.png",
       p_sp_darjeeling,
       width = 10, height = 8,
       units = "in", dpi = 300,bg = "white")
