## Provenance: Gemini (original), Claude (revision). See inst/PROVENANCE.md.
set.seed(1)
se <- makeMockSpatial()
analyzed <- suppressMessages(spatialLagGLM(se, ~condition, radius = 2.5))

test_that("plotSpatialLag returns ggplot objects for both metrics", {
    p1 <- plotSpatialLag(analyzed, "Gene1", "expression")
    p2 <- plotSpatialLag(analyzed, "Gene1", "spatial_lag")
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    ## The two metrics must show different values.
    expect_false(isTRUE(all.equal(p1$data$Value, p2$data$Value)))
})

test_that("plotSpatialLag facets by sample", {
    two <- suppressMessages(
        spatialLagGLM(makeTwoSampleMock(), ~1, radius = 2.5)
    )
    p <- plotSpatialLag(two, "Gene1", "spatial_lag")
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("plotSpatialLag validates its input", {
    expect_error(plotSpatialLag(1, "Gene1"), "SpatialExperiment")
    expect_error(plotSpatialLag(analyzed, "NoSuchGene"), "gene_id")
    expect_error(plotSpatialLag(analyzed, c("Gene1", "Gene2")), "gene_id")
    expect_error(plotSpatialLag(analyzed, "Gene1", "bogus"), "arg")
    expect_error(plotSpatialLag(analyzed, "Gene1", assay.type = "x"), "assay")
    expect_error(
        plotSpatialLag(se, "Gene1", "spatial_lag"), "Run spatialLagGLM"
    )
})

test_that("plotLagVolcano returns a ggplot and validates its input", {
    expect_s3_class(plotLagVolcano(analyzed), "ggplot")
    expect_error(plotLagVolcano(1), "SpatialExperiment")
    expect_error(plotLagVolcano(se), "Run spatialLagGLM")
    expect_error(plotLagVolcano(analyzed, alpha = 2), "alpha")
    expect_error(plotLagVolcano(analyzed, effect_cutoff = -1), "effect_cutoff")
})

test_that("runSpatialLagExplorer validates its input", {
    expect_error(runSpatialLagExplorer(1), "SpatialExperiment")
})
