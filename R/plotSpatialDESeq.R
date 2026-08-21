#' Plot Spatial Differential Expression Results
#'
#' @param se An analyzed SpatialExperiment object.
#' @param gene_id Character. Gene target to visualize.
#' @param metric Character. Choose "expression" or "spatial_effect".
#' @return A ggplot2 object.
#' @export
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay rowData
#' @import ggplot2
#' @export
#' @examples
#' se <- make_mock_spatial()
#' analyzed <- spatialDESeq(se, design = ~ condition, radius = 2.5)
#' plotSpatialDESeq(analyzed, gene_id = "Gene1", metric = "expression")
plotSpatialDESeq <- function(se, gene_id, metric = c("expression", "spatial_effect")) {
    if (!is(se, "SpatialExperiment")) stop("Input must be a SpatialExperiment object.")
    metric <- match.arg(metric)
    if (!gene_id %in% rownames(se)) stop("Gene not found in object.")
    
    coords <- as.data.frame(SpatialExperiment::spatialCoords(se))
    colnames(coords) <- c("X", "Y")
    gene_counts <- SummarizedExperiment::assay(se, "counts")[gene_id, ]
    plot_data <- cbind(coords, Expression = gene_counts)
    
    if (metric == "expression") {
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = X, y = Y, color = Expression)) +
            ggplot2::geom_point(size = 2, alpha = 0.8) +
            ggplot2::scale_color_viridis_c(option = "magma") +
            ggplot2::labs(title = paste("Spatial Expression:", gene_id), color = "Counts")
    } else {
        r_data <- as.data.frame(SummarizedExperiment::rowData(se))
        if (!"spatial_effect_size" %in% colnames(r_data)) stop("Run spatialDESeq first.")
        
        effect_size <- r_data[gene_id, "spatial_effect_size"]
        p_adj <- r_data[gene_id, "adjusted_p_value"]
        
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = X, y = Y, color = Expression)) +
            ggplot2::geom_point(size = 2, alpha = 0.8) +
            ggplot2::scale_color_viridis_c(option = "mako") +
            ggplot2::labs(title = paste("Spatial Model Context:", gene_id),
                          subtitle = sprintf("Effect: %.3f | Adj. p-val: %.3e", effect_size, p_adj),
                          color = "Counts")
    }
    
    p <- p + ggplot2::theme_minimal() + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.text = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank()) + 
        ggplot2::coord_fixed()
    return(p)
}
