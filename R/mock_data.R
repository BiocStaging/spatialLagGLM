## Provenance: Gemini (original), Claude (revision). See inst/PROVENANCE.md.
#' Generate Mock Spatial Transcriptomics Data
#'
#' Simulates a 30 x 30 grid of spots and 50 negative binomial genes. Spots
#' within 7 units of the grid centre are labelled \code{"TumorCore"}, the
#' rest \code{"Margin"}, and \code{Gene1} is up-regulated in the core.
#'
#' @return A \code{\link[SpatialExperiment]{SpatialExperiment}} object with
#'     synthetic counts and coordinates.
#' @export
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom stats rnbinom
#' @examples
#' se <- makeMockSpatial()
makeMockSpatial <- function() {
    grid_edge <- 30L
    coords <- expand.grid(px = seq_len(grid_edge), py = seq_len(grid_edge))
    num_spots <- nrow(coords)
    num_genes <- 50L

    center <- (grid_edge) / 2
    dist_to_center <- sqrt((coords$px - center)^2 + (coords$py - center)^2)
    condition <- ifelse(dist_to_center <= 7, "TumorCore", "Margin")

    counts <- matrix(
        stats::rnbinom(num_genes * num_spots, mu = 10, size = 2),
        nrow = num_genes, ncol = num_spots
    )
    core_spots <- which(condition == "TumorCore")
    counts[1, core_spots] <- stats::rnbinom(
        length(core_spots),
        mu = 120, size = 5
    )

    rownames(counts) <- paste0("Gene", seq_len(num_genes))
    colnames(counts) <- paste0("Spot", seq_len(num_spots))

    SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        spatialCoords = as.matrix(coords),
        colData = data.frame(
            condition = factor(condition),
            sample_id = "sample_01"
        )
    )
}
