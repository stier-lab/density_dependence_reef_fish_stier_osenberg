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

This repository contains the full data synthesis and analysis pipeline for a meta-analysis of density-dependent mortality in reef fishes. Using digitized data from >30 substudies across tropical coral reefs, we fit nonlinear mortality models to estimate density-independent (α) and density-dependent (β) mortality, examined ecological and phylogenetic drivers of β, and evaluated methodological effects.

This repository supports the manuscript:

**“Widespread Heterogeneity in Density-Dependent Mortality of Nearshore Fishes”**  
Stier & Osenberg, *in review*

---

## Repository Structure

```
code/           # All R scripts to load, analyze, and model the data
data/           # Raw digitized data, trait and covariate tables, phylogenies, metadata
output/         # Final model outputs and compiled estimates of α and β
results/        # Intermediate results, diagnostics, model comparisons
figures/        # Final visualizations used in the manuscript
```

---

## Quick Start

1. Clone the repository:

```bash
git clone https://github.com/[your-username]/density_dependence_reef_fish_stier_osenberg.git
cd density_dependence_reef_fish_stier_osenberg
```

2. Open R and set the working directory to the repository root.

3. Run the master script:

```r
source("code/00_run_all.R")
```

4. Required R packages are listed and loaded in `code/0_libraries.R`.

---

## Code Overview

| Script | Description |
|--------|-------------|
| `code/0_libraries.R` | Loads required R packages |
| `code/00_run_all.R` | Master script to run the entire pipeline |
| `code/1_data_phylogeny_loading.R` | Loads raw datasets and phylogenetic trees |
| `code/2_beta_estimate.R` | Fits nonlinear models to estimate α and β |
| `code/3_phylogenetic_analysis.R` | Tests for phylogenetic signal in β |
| `code/4_predators_pairedpredators.R` | Compares β with/without predators |
| `code/5_model_selection_comparison.R` | Trait-based model selection and meta-regression |
| `code/6_bivariate_plots_predictors.R` | Generates bivariate β vs. predictor plots |

---

## Data Contents

### `data/` folder includes:
- `all_studies_looped-2024-09-11.csv`: Raw mortality observations extracted from published figures/tables
- `combined_results_2024-09-17.csv`: Final α and β estimates for each substudy
- `covariates-2024-09-30.csv`: Ecological and methodological traits per substudy
- `unique_species_studies.xlsx`: Life history traits per species × substudy
- `manual_densities.csv`: Manually extracted density values for digitization QA/QC
- `*.tre`: Phylogenetic trees from Rabosky et al. (2018), Siqueira et al. (2020), and OpenTree

All data files are documented with accompanying `.txt` metadata files (Dryad-compatible).

---

## Output

### `output/`
- Model-fitted α and β with variances
- Tables of covariate effects
- Species-level and study-level summaries

### `figures/`
- Final figures used in the manuscript

### `results/`
- AIC tables
- Phylogenetic signal results
- Model comparison diagnostics

---

## Data Sources and Citations

Primary data were extracted from published ecological studies (see manuscript Table S1).  
Phylogenies sourced from:

- Rabosky et al. (2018): [https://doi.org/10.1038/s41586-018-0273-1](https://doi.org/10.1038/s41586-018-0273-1)
- Siqueira et al. (2020): [https://doi.org/10.1038/s41467-020-16498-w](https://doi.org/10.1038/s41467-020-16498-w)
- Open Tree of Life: [https://tree.opentreeoflife.org](https://tree.opentreeoflife.org)

---

## License

- **Data**: [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/)
- **Code**: MIT License

---

## Citation

When using this repository or dataset, please cite:

> Stier, A.C., & Osenberg, C.W. (in review). Widespread heterogeneity in density-dependent mortality of nearshore fishes. *Ecology Letters*.

> NSF Grant OCE-1851510
