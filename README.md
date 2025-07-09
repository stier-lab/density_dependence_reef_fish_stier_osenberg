![Coral reef fish](images/AdobeStock_466179536.jpeg)

# Widespread Heterogeneity in Density-Dependent Mortality of Nearshore Fishes

**Authors**  
Adrian C. Stier (University of California, Santa Barbara)  
Craig W. Osenberg (University of Georgia)  

**Contact**  
Adrian Stier — astier@ucsb.edu  
Craig Osenberg — osenberg@uga.edu  

**Funding**  
National Science Foundation (NSF) Grant [OCE-1851510](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1851510)

---

## Overview

This repository contains the complete data synthesis and analysis pipeline for a global meta-analysis on density-dependent mortality in reef fishes. We estimated mortality parameters (α and β) from >30 ecological studies and explored how ecological traits, experimental methods, and phylogenetic history explain variation in density dependence.

This repository supports the manuscript:

> **Stier & Osenberg**  
> *Widespread Heterogeneity in Density-Dependent Mortality of Nearshore Fishes*  
> *Ecology Letters* (in review)

---

## Repository Structure

```
code/           # All R scripts for data loading, modeling, and figure generation
data/           # Raw data, phylogenies, covariates, and metadata
output/         # Model outputs, effect sizes, and summary tables
results/        # Intermediate results, AIC tables, phylogenetic diagnostics
figures/        # Final figures used in the manuscript
images/         # Supplemental images for README (e.g., reef fish photo)
```

---

## Quick Start

1. Clone the repository:

```bash
git clone https://github.com/[your-username]/density_dependence_reef_fish_stier_osenberg.git
cd density_dependence_reef_fish_stier_osenberg
```

2. Open R or RStudio and run:

```r
source("code/00_run_all.R")
```

3. All outputs will be saved to the `output/`, `results/`, and `figures/` directories.

---

## Code Overview

| Script                                 | Description                                                         |
|----------------------------------------|---------------------------------------------------------------------|
| `code/0_libraries.R`                   | Loads all required R packages                                       |
| `code/00_run_all.R`                    | Master script to reproduce full pipeline                            |
| `code/1_data_phylogeny_loading.R`      | Loads raw data and phylogenies                                      |
| `code/2_beta_estimate.R`               | Estimates α and β from nonlinear fits                               |
| `code/3_phylogenetic_analysis.R`       | Quantifies phylogenetic signal in β                                 |
| `code/4_predators_pairedpredators.R`   | Tests for predator effects on β                                     |
| `code/5_model_selection_comparison.R`  | Compares covariate-based and random-effects models                  |
| `code/6_bivariate_plots_predictors.R`  | Generates bivariate plots of predictors vs. β                       |

---

## Data Contents

All data files are in the `data/` directory and are accompanied by `.txt` metadata files in Dryad-compatible format. Key files include:

- `all_studies_looped-2024-09-11.csv`: Raw digitized mortality data
- `combined_results_2024-09-17.csv`: Final α and β estimates across substudies
- `covariates-2024-09-30.csv`: Study-level traits and methodological variables
- `manual_densities.csv`: Manually extracted densities used for QA
- `unique_species_studies.xlsx`: Species traits per study
- `1.newick.tre`: Phylogeny from Open Tree of Life
- `actinopt_12k_treePL.tre`: Time-calibrated tree from Rabosky et al. (2018)
- `Reef_fish_all.tacted.newick.tre`: Species-level reef fish tree from Siqueira et al. (2020)

Each file has an associated metadata file (e.g., `*_metadata.txt`) describing variables, units, sources, and context.

---

## Output Contents

### `output/`
- Final α and β estimates with variances
- Covariate model coefficients and confidence intervals
- Exportable tables used in manuscript

### `results/`
- Phylogenetic signal estimates (Pagel’s λ)
- Model selection (AICc tables)
- Model diagnostics and fits

### `figures/`
- Final figures (e.g., trait-β plots, phylogenetic tree overlays)

---

## Dependencies

This project was developed using R (≥ 4.2.0). Core packages include:

- `metafor` – meta-analysis models
- `ape`, `phytools` – phylogenetic manipulation and plotting
- `dplyr`, `tibble`, `tidyr` – data wrangling
- `ggplot2`, `patchwork` – plotting
- `nlme`, `MCMCglmm` – mixed and phylogenetic models
- `readr`, `readxl` – file I/O
- `here`, `glue`, `stringr` – utilities

To install all packages, run:

```r
source("code/0_libraries.R")
```

---

## Workflow Instructions

To fully reproduce the analysis:

```r
# 1. Clone and enter repo
git clone https://github.com/[your-username]/density_dependence_reef_fish_stier_osenberg.git
cd density_dependence_reef_fish_stier_osenberg

# 2. Open R and run:
source("code/00_run_all.R")
```

This will:

1. Load and clean data  
2. Estimate α and β  
3. Fit phylogenetic and trait-based models  
4. Output diagnostics, figures, and tables

---

## Data Provenance

- **Mortality data**: Digitized from published ecological studies (see manuscript Table S1)
- **Phylogenies**:
  - Rabosky et al. (2018): [https://doi.org/10.1038/s41586-018-0273-1](https://doi.org/10.1038/s41586-018-0273-1)
  - Siqueira et al. (2020): [https://doi.org/10.1038/s41467-020-16498-w](https://doi.org/10.1038/s41467-020-16498-w)
  - OpenTree: [https://tree.opentreeoflife.org/curator/study/view/ot_1592](https://tree.opentreeoflife.org/curator/study/view/ot_1592)

---

## License

- **Data**: [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/)
- **Code**: [MIT License](https://opensource.org/licenses/MIT)

---

## Citation

If using this dataset or codebase, please cite:

> Stier, A. C. & Osenberg, C. W. *Widespread heterogeneity in density-dependent mortality of nearshore fishes.* Ecology Letters (in review).

Upon publication, citation DOIs (manuscript + Dryad) will be added.

---

## ORCID IDs

- Adrian C. Stier — [0000-0002-4704-4145](https://orcid.org/0000-0002-4704-4145)

