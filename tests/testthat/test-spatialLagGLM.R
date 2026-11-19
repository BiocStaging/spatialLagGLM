set.seed(1)
se <- makeMockSpatial()

test_that("results are returned in rowData", {
    out <- suppressMessages(spatialLagGLM(se, ~condition, radius = 2.5))
    expect_s4_class(out, "SpatialExperiment")
    expect_equal(dim(out), dim(se))
    rd <- SummarizedExperiment::rowData(out)
    expect_true(all(
        c("spatial_effect_size", "p_value", "adjusted_p_value") %in%
            colnames(rd)
    ))
    expect_equal(S4Vectors::metadata(out)$spatialLagGLM$radius, 2.5)
    expect_true(all(rd$adjusted_p_value >= rd$p_value, na.rm = TRUE))
})

test_that("Gene1 is detected when the design does not explain its pattern", {
    out <- suppressMessages(spatialLagGLM(se, ~1, radius = 2.5))
    rd <- SummarizedExperiment::rowData(out)
    expect_lt(rd["Gene1", "adjusted_p_value"], 0.001)
    expect_gt(rd["Gene1", "spatial_effect_size"], 0)
    expect_equal(which.min(rd$p_value), 1L)
})

test_that("multi-sample objects run", {
    out <- suppressMessages(
        spatialLagGLM(makeTwoSampleMock(), ~1, radius = 2.5)
    )
    expect_false(anyNA(SummarizedExperiment::rowData(out)$p_value))
})

test_that("assay.type selects the assay", {
    se2 <- se
    SummarizedExperiment::assay(se2, "other") <-
        SummarizedExperiment::assay(se2, "counts")
    a <- suppressMessages(spatialLagGLM(se2, ~1, assay.type = "other"))
    b <- suppressMessages(spatialLagGLM(se2, ~1, assay.type = "counts"))
    expect_equal(
        SummarizedExperiment::rowData(a)$p_value,
        SummarizedExperiment::rowData(b)$p_value
    )
})

test_that("size.factors can be supplied and are validated", {
    lib <- colSums(SummarizedExperiment::assay(se, "counts"))
    sf <- lib / exp(mean(log(lib)))
    a <- suppressMessages(spatialLagGLM(se, ~1, size.factors = sf))
    b <- suppressMessages(spatialLagGLM(se, ~1))
    expect_equal(
        SummarizedExperiment::rowData(a)$p_value,
        SummarizedExperiment::rowData(b)$p_value
    )
    expect_error(spatialLagGLM(se, ~1, size.factors = sf[-1]), "size.factors")
    expect_error(spatialLagGLM(se, ~1, size.factors = -sf), "size.factors")
})

test_that("inputs are validated", {
    expect_error(spatialLagGLM(1, ~condition), "SpatialExperiment")
    expect_error(spatialLagGLM(se, "condition"), "formula")
    expect_error(spatialLagGLM(se, ~condition, radius = -1), "radius")
    expect_error(spatialLagGLM(se, ~condition, radius = 0), "radius")
    expect_error(spatialLagGLM(se, ~condition, radius = c(1, 2)), "radius")
    expect_error(spatialLagGLM(se, ~condition, radius = "a"), "radius")
    expect_error(spatialLagGLM(se, ~condition, radius = NA_real_), "radius")
    expect_error(
        spatialLagGLM(se, ~condition, assay.type = "nope"), "not found"
    )
    expect_error(spatialLagGLM(se, ~condition, assay.type = 1), "assay.type")
})

test_that("spots with zero counts are rejected", {
    se0 <- se
    SummarizedExperiment::assay(se0, "counts")[, 1] <- 0L
    expect_error(spatialLagGLM(se0, ~1), "zero total counts")
})

test_that("missing counts are rejected", {
    se_na <- se
    SummarizedExperiment::assay(se_na, "counts")[2, 1] <- NA_integer_
    expect_error(spatialLagGLM(se_na, ~1), "missing counts")
})

test_that("failed fits are reported", {
    ## A negative count makes the quasi-Poisson fit of that gene fail.
    se_bad <- se
    SummarizedExperiment::assay(se_bad, "counts")[2, 1] <- -1L
    expect_warning(
        out <- suppressMessages(spatialLagGLM(se_bad, ~1)),
        "gene fits failed"
    )
    expect_true(is.na(SummarizedExperiment::rowData(out)["Gene2", "p_value"]))
})
