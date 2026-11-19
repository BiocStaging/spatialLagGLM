#' Plot a Gene on the Tissue Coordinates
#'
#' Plots either the raw counts of a gene or its spatial lag (the weighted mean
#' of its neighbours' normalised counts, i.e. the covariate used by
#' \code{\link{spatialLagGLM}}). With several samples, each is drawn in its
#' own panel.
#'
#' @param se A \code{SpatialExperiment} analysed with
#'     \code{\link{spatialLagGLM}} (needed only for \code{"spatial_lag"}).
#' @param gene_id A single gene name.
#' @param metric \code{"expression"} for the raw counts, or
#'     \code{"spatial_lag"} for the spatial lag.
#' @param assay.type A string naming the assay with raw counts.
#' @return A ggplot2 object.
#' @export
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay assayNames colData rowData
#' @importFrom S4Vectors metadata
#' @examples
#' se <- makeMockSpatial()
#' analyzed <- spatialLagGLM(se, design = ~condition, radius = 2.5)
#' plotSpatialLag(analyzed, gene_id = "Gene1", metric = "expression")
#' plotSpatialLag(analyzed, gene_id = "Gene1", metric = "spatial_lag")
plotSpatialLag <- function(
  se, gene_id,
  metric = c("expression", "spatial_lag"), assay.type = "counts"
) {
    .checkSpatialExperiment(se)
    metric <- match.arg(metric)
    .checkAssayType(se, assay.type)
    .checkGeneId(se, gene_id)

    counts <- as.matrix(SummarizedExperiment::assay(se, assay.type))
    coord_mat <- SpatialExperiment::spatialCoords(se)
    sample_ids <- .sampleIds(se)
    coords <- as.data.frame(coord_mat[, seq_len(2L), drop = FALSE])
    colnames(coords) <- c("X", "Y")
    plot_data <- cbind(coords, sample_id = sample_ids)

    if (metric == "expression") {
        plot_data$Value <- counts[gene_id, ]
        p <- ggplot2::ggplot(
            plot_data,
            ggplot2::aes(x = X, y = Y, color = Value)
        ) +
            ggplot2::geom_point(size = 2, alpha = 0.8) +
            ggplot2::scale_color_viridis_c(option = "magma") +
            ggplot2::labs(
                title = paste("Spatial expression:", gene_id),
                color = "Counts"
            )
    } else {
        info <- S4Vectors::metadata(se)$spatialLagGLM
        row_data <- as.data.frame(SummarizedExperiment::rowData(se))
        if (is.null(info) ||
            !all(c("spatial_effect_size", "adjusted_p_value") %in%
                colnames(row_data))) {
            stop("Run spatialLagGLM() first.", call. = FALSE)
        }
        weights <- .buildSpatialWeights(coord_mat, sample_ids, info$radius)
        size_factors <- info$size.factors
        if (is.null(size_factors)) {
            size_factors <- .sizeFactors(counts)
        }
        norm_counts <- counts[gene_id, ] / size_factors
        plot_data$Value <- as.numeric(weights %*% norm_counts)

        p <- ggplot2::ggplot(
            plot_data,
            ggplot2::aes(x = X, y = Y, color = Value)
        ) +
            ggplot2::geom_point(size = 2, alpha = 0.8) +
            ggplot2::scale_color_viridis_c(option = "mako") +
            ggplot2::labs(
                title = paste("Spatial lag:", gene_id),
                subtitle = sprintf(
                    "Effect: %.3f | Adj. p-value: %.3e",
                    row_data[gene_id, "spatial_effect_size"],
                    row_data[gene_id, "adjusted_p_value"]
                ),
                color = "Lag"
            )
    }

    if (length(unique(sample_ids)) > 1L) {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(sample_id))
    }
    p + ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank()
        ) +
        ggplot2::coord_fixed()
}
