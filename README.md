# spatialLagGLM

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`spatialLagGLM` tests, gene by gene, whether the expression at a spot is
associated with the expression at neighbouring spots in spatial
transcriptomics data. It works on `SpatialExperiment` objects and stores its
results in `rowData()`.

For each gene it fits a quasi-Poisson regression with a library-size offset,
optional design terms from `colData`, and the spatial lag of the gene's own
normalised counts as an extra covariate. The lag is an inverse-distance
weighted mean over spots within a given radius, computed separately for each
sample. Weights are computed in C++ (RcppArmadillo) and genes are fitted in
parallel with BiocParallel.

**Note on the statistics.** The null hypothesis is that the lag coefficient is
zero. The model is not a spatial autoregressive model, and the p-values are
exploratory: in simulation they were anti-conservative. See the vignette
sections on calibration and limitations before interpreting results. For a
well-established spatially variable gene test, consider SPARK-X, nnSVG or
SpatialDE.

## Installation

Once the package is accepted by Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("spatialLagGLM")
```

The development version from GitHub:

```r
BiocManager::install("sahalpaladan/spatialLagGLM")
```

## Quick start

```r
library(spatialLagGLM)

# Simulated 30 x 30 grid, 50 genes; Gene1 is high in a central region
se <- makeMockSpatial()

# Test the spatial lag, adjusting for `condition`
res <- spatialLagGLM(se, design = ~condition, radius = 2.5)
head(SummarizedExperiment::rowData(res))

# Plot a gene and the summary of all genes
plotSpatialLag(res, gene_id = "Gene1", metric = "spatial_lag")
plotLagVolcano(res)

# Interactive explorer (requires the shiny package)
runSpatialLagExplorer(res)
```

`Gene1` is explained by `condition` in this example, so it is not reported
with `design = ~condition`; use `design = ~1` to see it. The vignette explains
this and how to choose `radius` and `design`.

## Real data

A worked example on 10x Visium human DLPFC data from `STexampleData`, including
how to supply size factors when analysing a subset of genes, is in the
vignette:

```r
vignette("spatialLagGLM", package = "spatialLagGLM")
```

## Reporting problems

Please open an issue at
<https://github.com/sahalpaladan/spatialLagGLM/issues>.

## License

MIT. See the `LICENSE` file.
