## Two interleaved samples ("A" = odd spots, "B" = even spots), so a mix-up
## between samples would be detectable.
makeTwoSampleMock <- function() {
    se <- makeMockSpatial()
    SpatialExperiment::SpatialExperiment(
        assays = list(counts = SummarizedExperiment::assay(se, "counts")),
        spatialCoords = SpatialExperiment::spatialCoords(se),
        colData = S4Vectors::DataFrame(
            condition = se$condition,
            sample_id = rep(c("A", "B"), length.out = ncol(se))
        )
    )
}
