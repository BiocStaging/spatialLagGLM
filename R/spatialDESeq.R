#' Multithreaded Spatial Differential Expression Analysis
#'
#' @param se A SpatialExperiment object.
#' @param design A formula specifying experimental variables (e.g., ~ condition).
#' @param radius Numeric. Distance radius for local spatial neighbors.
#' @param BPPARAM A BiocParallelParam object defining core backend allocation.
#' @export
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom SummarizedExperiment assay rowData rowData<- colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom stats model.matrix glm p.adjust quasipoisson
#' @import methods
#' @importFrom Rcpp evalCpp
#' @export
#' @useDynLib spatialDESeq, .registration = TRUE
#' @examples
#' se <- make_mock_spatial()
#' analyzed <- spatialDESeq(se, design = ~ condition, radius = 2.5)
#' @return A \code{SpatialExperiment} object with analysis results appended to \code{rowData}.
spatialDESeq <- function(se, design, radius = 2.5, BPPARAM = BiocParallel::bpparam()) {
    
    if (!is(se, "SpatialExperiment")) stop("Input must be a SpatialExperiment object.")
    
    counts <- as.matrix(SummarizedExperiment::assay(se, "counts"))
    coords <- SpatialExperiment::spatialCoords(se)
    metadata <- as.data.frame(SummarizedExperiment::colData(se))
    design_matrix <- stats::model.matrix(design, data = metadata)
    
    message("--> Calculating library-size normalization factors...")
    lib_sizes <- colSums(counts)
    size_factors <- lib_sizes / exp(mean(log(lib_sizes[lib_sizes > 0])))
    
    message("--> Computing C++ physical distance matrices...")
    weight_matrix <- rcpp_spatial_weights(coords, radius)
    
    message("--> Launching parallel statistical engine across cores...")
    num_genes <- nrow(counts)
    
    results <- BiocParallel::bplapply(seq_len(num_genes), function(i) {
        gene_counts <- counts[i, ]
        norm_counts <- gene_counts / size_factors
        spatial_covariate <- as.vector(weight_matrix %*% norm_counts)
        
        full_data <- data.frame(counts = gene_counts, design_matrix, spatial_cov = spatial_covariate)
        
        fit <- tryCatch({
            stats::glm(counts ~ . - counts + offset(log(size_factors)), 
                       data = full_data, family = stats::quasipoisson())
        }, error = function(e) { NULL })
        
        if (is.null(fit)) return(c(p_val = NA, estimate = NA))
        
        coef_summary <- summary(fit)$coefficients
        if ("spatial_cov" %in% rownames(coef_summary)) {
            return(c(p_val = coef_summary["spatial_cov", "Pr(>|t|)"], 
                     estimate = coef_summary["spatial_cov", "Estimate"]))
        } else {
            return(c(p_val = NA, estimate = NA))
        }
    }, BPPARAM = BPPARAM)
    
    res_df <- do.call(rbind, results)
    SummarizedExperiment::rowData(se)$spatial_effect_size <- res_df[, "estimate"]
    SummarizedExperiment::rowData(se)$p_value <- res_df[, "p_val"]
    SummarizedExperiment::rowData(se)$adjusted_p_value <- stats::p.adjust(res_df[, "p_val"], method = "BH")
    
    message("--> Done! Spatial evaluations appended to object metadata.")
    return(se)
}
