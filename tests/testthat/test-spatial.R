library(testthat)
library(spatialDESeq)

test_that("Full structural package calculations execute perfectly", {
    se <- make_mock_spatial()
    out <- spatialDESeq(se, design = ~ condition, radius = 2.0)
    expect_s4_class(out, "SpatialExperiment")
    expect_true("adjusted_p_value" %in% colnames(SummarizedExperiment::rowData(out)))
})
