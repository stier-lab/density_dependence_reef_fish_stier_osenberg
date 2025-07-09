###############################################################################
# File: 1_data_phylogeny_loading.R
# Purpose: Load and preprocess trait metadata and phylogenetic trees for analysis
# Author: Adrian Stier
# Date: 2025-07-08
###############################################################################

# -----------------------------------------------------------------------------
# 1. Dependencies & Preliminary Sources
# -----------------------------------------------------------------------------
source(here::here("code", "0_libraries.R"))                       # Load required packages
source(here::here("code", "4_combine_params_covariates.R"))     # all_dat2 object preparation
source(here::here("code", "5_transformations.R"))               # Data transformations

# -----------------------------------------------------------------------------
# 2. Phylogenetic Trees Loading
# -----------------------------------------------------------------------------
reef_fish_tree   <- ape::read.tree(here::here("data", "Reef_fish_all.tacted.newick.tre"))
actino_tree      <- ape::read.tree(here::here("data", "actinopt_12k_treePL.tre"))
custom_tree_12k  <- ape::read.tree(here::here("data", "1.newick.tre"))
fish_phylo_tree  <- fishtree_phylogeny()  # Latest fish phylogeny from fishtree package

# -----------------------------------------------------------------------------
# 3. Load and Merge Species Metadata
# -----------------------------------------------------------------------------
meta_phy_df <- readxl::read_excel(
  here::here("data", "unique_species_studies.xlsx")
) %>%
  dplyr::rename(
    max_length_cm = `max_len (cm)`,
    max_weight_g  = `Max_wt (ga)`
  )

all_dat2 <- all_dat2 %>%
  dplyr::left_join(
    meta_phy_df,
    by = c("g_sp", "study_num", "substudy_num")
  ) %>%
  dplyr::mutate(
    max_length_density = max_length_cm * mean_density
  )

# -----------------------------------------------------------------------------
# 4. Unique Species-Study Combinations
# -----------------------------------------------------------------------------
species_study_df <- all_dat2 %>%
  dplyr::distinct(g_sp, study_num, substudy_num) %>%
  dplyr::arrange(g_sp, study_num, substudy_num)

readr::write_csv(
  species_study_df,
  here::here("data", "unique_species_studies.csv")
)
message("Saved unique species-study combinations: ", nrow(species_study_df), " rows")

# -----------------------------------------------------------------------------
# 5. Summary Statistics
# -----------------------------------------------------------------------------
n_unique_species <- all_dat2 %>%
  dplyr::pull(g_sp) %>%
  unique() %>%
  length()
message("Number of unique species: ", n_unique_species)

n_studies    <- dplyr::n_distinct(all_dat2$study_num)
n_substudies <- dplyr::n_distinct(all_dat2$substudy_num)

message("Number of unique studies:   ", n_studies)
message("Number of unique substudies:", n_substudies)

# -----------------------------------------------------------------------------
# 6. Publication Information Tables
# -----------------------------------------------------------------------------
study_substudy_info_df <- all_dat2 %>%
  dplyr::select(
    study_num,
    substudy_num,
    Authors,
    Article.Title,
    Source.Title,
    Publication.Year
  ) %>%
  dplyr::distinct()

study_info_df <- all_dat2 %>%
  dplyr::select(
    study_num,
    Authors,
    Article.Title,
    Source.Title,
    Publication.Year
  ) %>%
  dplyr::distinct()

readr::write_csv(
  study_substudy_info_df,
  here::here("output", "unique_study_substudy_info.csv")
)
readr::write_csv(
  study_info_df,
  here::here("output", "unique_study_info.csv")
)

# Optional: Render a formatted table for reporting
study_substudy_info_df %>%
  gt::gt() %>%
  gt::tab_header(
    title    = "Study / Substudy → Publication Details",
    subtitle = "Authors • Article Title • Journal • Year"
  ) %>%
  gt::cols_label(
    study_num        = "Study #",
    substudy_num     = "Substudy #",
    Authors          = "Authors",
    Article.Title    = "Article Title",
    Source.Title     = "Journal",
    Publication.Year = "Year"
  ) -> study_table

gtsave(
  study_table,
  filename = here::here("output", "study_substudy_info.html")
)

# -----------------------------------------------------------------------------
# 7. Meta-Analysis Citation Table
# -----------------------------------------------------------------------------
# Dependencies: dplyr, stringr, gt (ensure gt is loaded)
library(dplyr)
library(stringr)
library(gt)

# Helper: Standardize author names to "Surname, INITIALS"
fix_author_case <- function(name) {
  name <- str_trim(name)
  parts <- str_split(name, ",\\s*")[[1]]
  surname <- str_to_title(parts[1])
  initials <- if (length(parts) > 1) {
    str_replace_all(parts[2], "\\b\\w+\\b", function(x) toupper(x))
  } else {
    ""
  }
  str_trim(paste0(surname, if (initials != "") paste0(", ", initials)))
}

# Helper: Ensure proper capitalization in author names
clean_author_name <- function(name) {
  # Trim and split by spaces or commas
  parts <- unlist(str_split(str_trim(name), "\\s+|,\\s*| and "))
  cleaned <- sapply(parts, function(word) {
    if (str_detect(word, "^[A-Z]{2,}$")) word else str_to_sentence(word)
  })
  # Recombine preserving initials
  str_replace_all(name, "\\b[A-Za-z]+\\b", function(x) {
    if (str_detect(x, "^[A-Z]{2,}$")) x else str_to_sentence(x)
  })
}

# Build the meta-analysis summary table
meta_table <- all_dat2 %>%
  mutate(
    Author   = Authors,
    Year     = Publication.Year,
    Citation = ifelse(
      !is.na(Authors) & !is.na(`Article.Title`) & !is.na(`Source.Title`) & !is.na(Publication.Year),
      paste0(
        str_to_title(Authors), " (", Year, "). ",
        str_to_sentence(`Article.Title`), ". ",
        "*", str_to_title(`Source.Title`), "*."
      ),
      NA_character_
    )
  ) %>%
  select(
    Study        = study_num,
    Substudy     = substudy_num,
    Description  = data.source..or.reason.for.exclusion.,
    Beta         = betanls2_raw_cm,
    Variance     = betanlsvar_raw_cm,
    Author,
    Year,
    Citation,
    DOI,
    g_sp,
    PairedStudy  = paired_pred,
    ExptObsPair  = expt_obs_pairs,
    Predators    = predators,
    Duration     = duration,
    Exp_obs      = expt_obs,
    SizeStart    = size_start,
    MaxLen       = max_len_cm,
    mean_density,
    family       = family
  ) %>%
  distinct() %>%
  arrange(Year, Author) %>%
  mutate(
    Beta     = round(Beta, 3),
    Variance = round(Variance, 5)
  )

# Generate and save a formatted GT table
meta_table %>%
  gt() %>%
  tab_header(
    title    = md("**Summary of Studies Included in Meta-Analysis**"),
    subtitle = "Effect Size, Variance, and Metadata by Substudy"
  ) %>%
  fmt_number(columns = c("Beta", "Variance"), decimals = 3) %>%
  cols_label(
    Study       = "Study #",
    Substudy    = "Substudy #",
    Beta        = "Effect Size (β)",
    Variance    = "Variance (asinh)",
    Author      = "First Author(s)",
    Year        = "Year",
    Citation    = "Full Citation",
    PairedStudy = "Paired Predatory Treatment",
    ExptObsPair = "Exp/Obs Pair ID"
  ) %>%
  tab_options(
    table.font.size         = px(13),
    column_labels.font.weight = "bold",
    heading.align           = "center"
  ) -> meta_gt_table

gtsave(
  meta_gt_table,
  filename = here::here("output", "meta_analysis_table.html")
)

# Write raw meta_table for external use
readr::write_csv(
  meta_table,
  here::here("figures", "meta_analysis_table.csv"),
  na = ""
)

# End of 1_data_phylogeny_loading.R
