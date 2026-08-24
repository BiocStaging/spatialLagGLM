# spatialDESeq

<!-- Badges Grid Layout -->
[![BiocCheck Status](https://shields.io)](https://github.com)
[![R-CMD-check](https://shields.io)](https://github.com)
[![Ecosystem](https://shields.io)](https://bioconductor.org)
[![License: MIT](https://shields.io)](https://opensource.org)

**spatialDESeq** introduces an accelerated, parallelized framework optimized for spatial transcriptomics datasets. It extends standard negative binomial generalized linear models (GLMs) by incorporating physical tissue spot coordinates as distance-weighted local smoothing covariates, compiled natively via a C++ Armadillo matrix engine.

## 🌟 Key Functional Enhancements
* **`RcppArmadillo` Engine**: Avoids slow native loops by executing physical tissue distance mapping within a compiled vector space.
* **Continuous Neighborhood Smoothing**: Avoids hard category constraints by evaluating coordinate transcript densities.
* **Bioconductor Container Integration**: Operates natively on standard `SpatialExperiment` structures.
* **Point-and-Click Dashboard**: Built-in interactive workspace browser for local execution.

---

## 💾 System Installation

Install the package directly from your R workspace console using the `devtools` engine:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install core infrastructure environments
BiocManager::install(c("SpatialExperiment", "SummarizedExperiment", "BiocParallel"))

# Install spatialDESeq directly from source
devtools::install_github("sahalpaladan/spatialDESeq")
```

---

## 🧬 Quick-Start Code Workflow

Evaluate a complete modeling run using the built-in synthetic tissue generator matrix:

```r
library(spatialDESeq)
library(SpatialExperiment)
library(SummarizedExperiment)

# 1. Generate a mock tissue slide grid (50 genes, 900 spots)
se_raw <- make_mock_spatial()

# 2. Run the C++ accelerated spatial autoregressive solver
# Targets adjacent spatial coordinates within a radius footprint
se_analyzed <- spatialDESeq(se_raw, design = ~ condition, radius = 2.5)

# 3. Inspect calculated metrics appended to your object
head(rowData(se_analyzed))

# 4. Generate spatial context profile plots
plotSpatialDESeq(se_analyzed, gene_id = "Gene1", metric = "spatial_effect")

# 5. Open the interactive point-and-click visualization dashboard board
runExplorer(se_analyzed)
```

---

## 📊 Run Analytics on Real 10x Genomics Visium Data

```r
library(STexampleData)

# Load real public human brain tissue coordinates
se_real <- STexampleData::Visium_humanDLPFC()
rownames(se_real) <- make.names(rownames(se_real))

# Process spatial matrix engine (Radius tuned to large pixel grids)
se_results <- spatialDESeq(se_real[1:100, ], design = ~ ground_truth, radius = 150)

# Display multi-gene summary Volcano metrics
plotVolcano(se_results, alpha = 0.05, effect_cutoff = 0.5)
```

---

## 📋 Software Development Guidelines

### Bug Reports & Contributions
If you identify scaling bottlenecks, please submit an entry directly to the [GitHub Issues Portal](https://github.com/issues). 

### Licensing
This software ecosystem operates under the strict guidelines of the **MIT License**. See the `LICENSE` file for more details.
