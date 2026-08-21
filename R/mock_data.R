#' Generate Mock Spatial Transcriptomics Data
#' @return A SpatialExperiment object with synthetic counts and coordinates.
#' @export
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stats rnbinom
#' @examples
#' se <- make_mock_spatial()
make_mock_spatial <- function() {
    grid_edge <- 30
    coords <- expand.grid(x = seq_len(grid_edge), y = seq_len(grid_edge))
    colnames(coords) <- c("px", "py")
    num_spots <- nrow(coords)
    num_genes <- 50
    
    center_x <- 15
    center_y <- 15
    distances_to_center <- sqrt((coords$px - center_x)^2 + (coords$py - center_y)^2)
    condition <- ifelse(distances_to_center <= 7, "TumorCore", "Margin")
    
    counts <- matrix(stats::rnbinom(num_genes * num_spots, mu = 10, size = 2), 
                     nrow = num_genes, ncol = num_spots)
    
    core_spots <- which(condition == "TumorCore")
    counts[1, core_spots] <- stats::rnbinom(length(core_spots), mu = 120, size = 5) 
    
    rownames(counts) <- paste0("Gene", seq_len(num_genes))
    colnames(counts) <- paste0("Spot", seq_len(num_spots))
    
    se <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        spatialCoords = as.matrix(coords),
        colData = data.frame(condition = factor(condition), sample_id = "sample_01")
    )
    
    return(se)
}
