test_that("weights match a hand-computed example", {
    ## Points on a line at x = 0, 1, 2; all pairs lie within radius 2.5.
    ## The C++ uses weights 1/d, row-normalised, self excluded:
    ## row 1: d = 1, 2   -> 2/3, 1/3
    ## row 2: d = 1, 1   -> 1/2, 1/2
    ## row 3: d = 2, 1   -> 1/3, 2/3
    coords <- cbind(x = c(0, 1, 2), y = c(0, 0, 0))
    expected <- matrix(
        c(
            0, 2 / 3, 1 / 3,
            1 / 2, 0, 1 / 2,
            1 / 3, 2 / 3, 0
        ),
        nrow = 3, byrow = TRUE
    )
    w <- spatialLagGLM:::rcpp_spatial_weights(coords, 2.5)
    expect_equal(unname(as.matrix(w)), expected)
})

test_that("rows with neighbours sum to one and the diagonal is zero", {
    coords <- as.matrix(expand.grid(x = 1:5, y = 1:5))
    w <- as.matrix(spatialLagGLM:::rcpp_spatial_weights(coords, 1.5))
    expect_equal(unname(diag(w)), rep(0, nrow(coords)))
    expect_equal(unname(rowSums(w)), rep(1, nrow(coords)))
})

test_that("weights never link spots from different samples", {
    se <- makeTwoSampleMock()
    sid <- SummarizedExperiment::colData(se)$sample_id
    w <- spatialLagGLM:::.buildSpatialWeights(
        SpatialExperiment::spatialCoords(se), sid,
        radius = 2.5
    )
    expect_equal(dim(w), c(ncol(se), ncol(se)))
    nz <- Matrix::which(w != 0, arr.ind = TRUE)
    expect_gt(nrow(nz), 0)
    expect_true(all(sid[nz[, 1]] == sid[nz[, 2]]))
    expect_equal(unname(Matrix::rowSums(w)), rep(1, ncol(se)))
})

test_that("an isolated spot gets a zero row", {
    coords <- cbind(x = c(0, 1, 100), y = c(0, 0, 0))
    w <- spatialLagGLM:::.buildSpatialWeights(coords, rep("A", 3), 2)
    expect_equal(unname(Matrix::rowSums(w)), c(1, 1, 0))
})
