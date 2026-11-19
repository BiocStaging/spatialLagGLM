## Quiets R CMD check about non-standard evaluation in ggplot2 aesthetics.
utils::globalVariables(c(
    "X", "Y", "Value", "sample_id",
    "spatial_effect_size", "adjusted_p_value", "Significance"
))

.checkSpatialExperiment <- function(se) {
    if (!methods::is(se, "SpatialExperiment")) {
        stop("'se' must be a SpatialExperiment object.", call. = FALSE)
    }
    invisible(TRUE)
}

.checkAssayType <- function(se, assay.type) {
    if (!is.character(assay.type) || length(assay.type) != 1L ||
        is.na(assay.type)) {
        stop("'assay.type' must be a single string.", call. = FALSE)
    }
    if (!assay.type %in% SummarizedExperiment::assayNames(se)) {
        stop(sprintf("assay '%s' not found in 'se'.", assay.type),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.checkRadius <- function(radius) {
    if (!is.numeric(radius) || length(radius) != 1L || !is.finite(radius) ||
        radius <= 0) {
        stop("'radius' must be a single positive number.", call. = FALSE)
    }
    invisible(TRUE)
}

.checkGeneId <- function(se, gene_id) {
    if (!is.character(gene_id) || length(gene_id) != 1L || is.na(gene_id) ||
        !gene_id %in% rownames(se)) {
        stop("'gene_id' must be a single gene name present in 'se'.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

## Size factors: library size scaled by its geometric
## mean. Spots with no counts have no defined offset, so they are rejected.
.sizeFactors <- function(counts) {
    lib_sizes <- colSums(counts)
    if (anyNA(lib_sizes) || any(lib_sizes <= 0)) {
        stop("spots with zero total counts or missing counts found; ",
            "remove them first.",
            call. = FALSE
        )
    }
    lib_sizes / exp(mean(log(lib_sizes)))
}

.checkSizeFactors <- function(size.factors, n) {
    if (!is.numeric(size.factors) || length(size.factors) != n ||
        anyNA(size.factors) || any(!is.finite(size.factors)) ||
        any(size.factors <= 0)) {
        stop("'size.factors' must be NULL or a positive numeric vector ",
            "with one value per spot.",
            call. = FALSE
        )
    }
    unname(size.factors)
}

.sampleIds <- function(se) {
    sid <- SummarizedExperiment::colData(se)[["sample_id"]]
    if (is.null(sid)) {
        return(rep("sample01", ncol(se)))
    }
    as.character(sid)
}

## Sparse n x n weight matrix, in the column order of the input. Weights are
## computed separately for each sample, so spots from different sections are
## never neighbours.
.buildSpatialWeights <- function(coords, sample_ids, radius) {
    coords <- as.matrix(coords)
    storage.mode(coords) <- "double"
    ii <- list()
    jj <- list()
    xx <- list()
    for (s in unique(sample_ids)) {
        idx <- which(sample_ids == s)
        w <- rcpp_spatial_weights(coords[idx, , drop = FALSE], radius)
        w[!is.finite(w)] <- 0
        nz <- which(w != 0, arr.ind = TRUE)
        ii[[s]] <- idx[nz[, 1L]]
        jj[[s]] <- idx[nz[, 2L]]
        xx[[s]] <- w[nz]
    }
    Matrix::sparseMatrix(
        i = unlist(ii, use.names = FALSE),
        j = unlist(jj, use.names = FALSE),
        x = unlist(xx, use.names = FALSE),
        dims = c(nrow(coords), nrow(coords))
    )
}
