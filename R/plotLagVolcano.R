## Provenance: Gemini (original), Claude (revision). See inst/PROVENANCE.md.
#' Volcano Plot of Spatial Lag Results
#'
#' @param se A \code{SpatialExperiment} analysed with
#'     \code{\link{spatialLagGLM}}.
#' @param alpha A single number in (0, 1]: adjusted p-value cutoff.
#' @param effect_cutoff A single non-negative number: effect size cutoff.
#' @return A ggplot2 object.
#' @export
#' @importFrom SummarizedExperiment rowData
#' @examples
#' se <- makeMockSpatial()
#' analyzed <- spatialLagGLM(se, design = ~condition, radius = 2.5)
#' plotLagVolcano(analyzed)
plotLagVolcano <- function(se, alpha = 0.05, effect_cutoff = 0.5) {
    .checkSpatialExperiment(se)
    if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
        alpha <= 0 || alpha > 1) {
        stop("'alpha' must be a single number in (0, 1].", call. = FALSE)
    }
    if (!is.numeric(effect_cutoff) || length(effect_cutoff) != 1L ||
        is.na(effect_cutoff) || effect_cutoff < 0) {
        stop("'effect_cutoff' must be a single non-negative number.",
            call. = FALSE
        )
    }
    res_df <- as.data.frame(SummarizedExperiment::rowData(se))
    if (!all(c("spatial_effect_size", "adjusted_p_value") %in%
        colnames(res_df))) {
        stop("Run spatialLagGLM() first.", call. = FALSE)
    }

    keep <- !is.na(res_df$adjusted_p_value) &
        !is.na(res_df$spatial_effect_size)
    res_df <- res_df[keep, , drop = FALSE]
    sig <- res_df$adjusted_p_value < alpha
    res_df$Significance <- "Not Significant"
    res_df$Significance[sig & res_df$spatial_effect_size > effect_cutoff] <-
        "Spatially Enriched"
    res_df$Significance[sig & res_df$spatial_effect_size < -effect_cutoff] <-
        "Spatially Depleted"

    ggplot2::ggplot(
        res_df,
        ggplot2::aes(
            x = spatial_effect_size, y = -log10(adjusted_p_value),
            color = Significance
        )
    ) +
        ggplot2::geom_point(alpha = 0.7, size = 1.8) +
        ggplot2::scale_color_manual(values = c(
            "Not Significant" = "grey70",
            "Spatially Enriched" = "firebrick3",
            "Spatially Depleted" = "dodgerblue3"
        )) +
        ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
        ggplot2::geom_vline(
            xintercept = c(-effect_cutoff, effect_cutoff),
            linetype = "dashed"
        ) +
        ggplot2::labs(
            title = "Spatial lag effect volcano plot",
            x = "Spatial effect size", y = "-log10(adjusted p-value)"
        ) +
        ggplot2::theme_classic()
}
