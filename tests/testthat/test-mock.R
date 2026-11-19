test_that("makeMockSpatial gives the documented object", {
    se <- makeMockSpatial()
    expect_s4_class(se, "SpatialExperiment")
    expect_equal(dim(se), c(50L, 900L))
    expect_equal(rownames(se)[1], "Gene1")
    expect_setequal(levels(se$condition), c("Margin", "TumorCore"))
    expect_equal(ncol(SpatialExperiment::spatialCoords(se)), 2L)
})
