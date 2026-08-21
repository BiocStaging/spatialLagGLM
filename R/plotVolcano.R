#' Generate Volcano Plot for Spatial Differential Expression
#'
#' @param se An analyzed SpatialExperiment object.
#' @param alpha Numeric. Significance value cutoff.
#' @param effect_cutoff Numeric. Threshold for effect markings.
#' @return A ggplot2 object.
#' @export
#' @importFrom SummarizedExperiment rowData
#' @import ggplot2
#' @examples
#' se <- make_mock_spatial()
#' analyzed <- spatialDESeq(se, design = ~ condition, radius = 2.5)
#' plotVolcano(analyzed)
plotVolcano <- function(se, alpha = 0.05, effect_cutoff = 0.5) {
    if (!is(se, "SpatialExperiment")) stop("Input must be a SpatialExperiment object.")
    res_df <- as.data.frame(SummarizedExperiment::rowData(se))
    if (!all(c("spatial_effect_size", "adjusted_p_value") %in% colnames(res_df))) stop("Run spatialDESeq first.")
    
    res_df <- res_df[!is.na(res_df$adjusted_p_value) & !is.na(res_df$spatial_effect_size), ]
    res_df$Significance <- "Not Significant"
    res_df$Significance[res_df$adjusted_p_value < alpha & res_df$spatial_effect_size > effect_cutoff] <- "Spatially Enriched"
    res_df$Significance[res_df$adjusted_p_value < alpha & res_df$spatial_effect_size < -effect_cutoff] <- "Spatially Depleted"
    
    p <- ggplot2::ggplot(res_df, ggplot2::aes(x = spatial_effect_size, y = -log10(adjusted_p_value), color = Significance)) +
        ggplot2::geom_point(alpha = 0.7, size = 1.8) +
        ggplot2::scale_color_manual(values = c("Not Significant" = "grey70", "Spatially Enriched" = "firebrick3", "Spatially Depleted" = "dodgerblue3")) +
        ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
        ggplot2::geom_vline(xintercept = c(-effect_cutoff, effect_cutoff), linetype = "dashed") +
        ggplot2::labs(title = "Multi-Gene Target Volcano Summary", x = "Spatial Effect Size", y = "-log10(Adj p-value)") +
        ggplot2::theme_classic()
    return(p)
}
