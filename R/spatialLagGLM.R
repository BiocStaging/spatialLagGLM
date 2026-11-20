## Provenance: Gemini (original), Claude (revision). See inst/PROVENANCE.md.
#' Quasi-Poisson Regression with a Spatially Lagged Covariate
#'
#' For every gene, fits a quasi-Poisson regression with a library-size offset,
#' the terms of \code{design}, and the spatial lag of the gene's own
#' size-factor-normalised counts as an extra covariate. The coefficient of the
#' spatial lag, its p-value and its Benjamini-Hochberg adjusted p-value are
#' stored in \code{rowData(se)}.
#'
#' @details
#' For spot \eqn{i} and a given gene, the model is
#' \deqn{\log E[y_i] = \log(s_i) + X_i \beta + \gamma \, l_i,}
#' where \eqn{s_i} is the size factor (library size divided by its geometric
#' mean), \eqn{X_i} the row of the design matrix, and
#' \eqn{l_i = \sum_j w_{ij} y_j / s_j} the spatial lag. The weights
#' \eqn{w_{ij}} are inverse distances between spots within \code{radius},
#' row-normalised, and are computed separately for each \code{sample_id}.
#'
#' The null hypothesis tested is \eqn{\gamma = 0}: given the offset and the
#' design terms, a spot's expression is not associated with the weighted mean
#' of its neighbours' normalised expression.
#'
#' Because the lag is built from the same gene's counts, this is not a
#' spatial autoregressive model, and the quasi-Poisson p-values have not been
#' shown to be calibrated for this test. Treat them as exploratory and, where
#' it matters, check them against a permutation of the spot coordinates.
#'
#' Genes whose fit fails or does not converge get \code{NA} results, and a
#' warning reports how many.
#'
#' The results are stored in \code{rowData(se)} as \code{spatial_effect_size}
#' (\eqn{\hat\gamma}), \code{p_value} and \code{adjusted_p_value}. The
#' radius, assay and size factors used are stored in
#' \code{metadata(se)$spatialLagGLM}.
#'
#' @param se A \code{\link[SpatialExperiment]{SpatialExperiment}} object.
#' @param design A formula using variables in \code{colData(se)}, e.g.
#'     \code{~ condition}.
#' @param radius A positive number: the distance within which spots are
#'     neighbours, in the units of \code{spatialCoords(se)}. For several
#'     samples, all samples should share the same units.
#' @param assay.type A string naming the assay with raw counts.
#' @param size.factors \code{NULL} (default) to compute size factors from
#'     the library sizes of \code{assay.type}, or a positive numeric vector
#'     with one value per spot. Supply this when \code{se} has been
#'     subsetted to a few genes, so that normalisation still reflects all
#'     genes.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#' @return \code{se} with \code{spatial_effect_size}, \code{p_value} and
#'     \code{adjusted_p_value} added to \code{rowData}.
#' @export
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom SummarizedExperiment assay assayNames colData rowData
#'     rowData<-
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats model.matrix glm p.adjust quasipoisson
#' @examples
#' se <- makeMockSpatial()
#' analyzed <- spatialLagGLM(se, design = ~condition, radius = 2.5)
#' head(SummarizedExperiment::rowData(analyzed))
spatialLagGLM <- function(
  se, design, radius = 2.5, assay.type = "counts",
  size.factors = NULL, BPPARAM = BiocParallel::bpparam()
) {
    .checkSpatialExperiment(se)
    if (!inherits(design, "formula")) {
        stop("'design' must be a formula.", call. = FALSE)
    }
    .checkRadius(radius)
    .checkAssayType(se, assay.type)

    counts <- as.matrix(SummarizedExperiment::assay(se, assay.type))
    metadata_df <- as.data.frame(SummarizedExperiment::colData(se))
    design_matrix <- stats::model.matrix(design, data = metadata_df)
    if (nrow(design_matrix) != ncol(counts)) {
        stop("'design' variables contain missing values.", call. = FALSE)
    }

    size_factors <- if (is.null(size.factors)) {
        .sizeFactors(counts)
    } else {
        .checkSizeFactors(size.factors, ncol(counts))
    }

    message("--> Computing per-sample spatial weights...")
    weights <- .buildSpatialWeights(
        SpatialExperiment::spatialCoords(se), .sampleIds(se), radius
    )
    ## One sparse product for all genes: spots x genes.
    spatial_lag <- as.matrix(weights %*% (t(counts) / size_factors))

    message("--> Fitting gene-wise models...")
    results <- BiocParallel::bplapply(
        seq_len(nrow(counts)), .fitGene,
        counts = counts, spatial_lag = spatial_lag,
        design_matrix = design_matrix, log_sf = log(size_factors),
        BPPARAM = BPPARAM
    )
    res <- do.call(rbind, results)

    n_failed <- sum(is.na(res[, "p_value"]))
    if (n_failed > 0L) {
        warning(sprintf(
            "%d of %d gene fits failed or did not converge; results are NA.",
            n_failed, nrow(res)
        ), call. = FALSE)
    }

    row_data <- SummarizedExperiment::rowData(se)
    row_data[["spatial_effect_size"]] <- res[, "estimate"]
    row_data[["p_value"]] <- res[, "p_value"]
    row_data[["adjusted_p_value"]] <- stats::p.adjust(res[, "p_value"], "BH")
    SummarizedExperiment::rowData(se) <- row_data
    S4Vectors::metadata(se)$spatialLagGLM <- list(
        radius = radius, assay.type = assay.type,
        size.factors = size_factors
    )

    message("--> Done.")
    se
}

## Fits one gene. Returns c(estimate, p_value), both NA on failure.
.fitGene <- function(i, counts, spatial_lag, design_matrix, log_sf) {
    failed <- c(estimate = NA_real_, p_value = NA_real_)
    y <- counts[i, ]
    X <- cbind(design_matrix, spatial_cov = spatial_lag[, i])
    fit <- tryCatch(
        stats::glm(y ~ 0 + X, family = stats::quasipoisson(), offset = log_sf),
        error = function(e) NULL
    )
    if (is.null(fit) || !fit$converged) {
        return(failed)
    }
    coefs <- summary(fit)$coefficients
    if (!"Xspatial_cov" %in% rownames(coefs)) {
        return(failed)
    }
    c(
        estimate = coefs["Xspatial_cov", "Estimate"],
        p_value = coefs["Xspatial_cov", "Pr(>|t|)"]
    )
}
