# testing mapout
testdata <- matrix(c(rnorm(100*100, sd=2), rnorm(100*100, sd=5),
                  rnorm(100*100, sd=3), rnorm(100*100, sd=7)), nrow = 100)
rownames(testdata) <- as.character(1:nrow(testdata))
obj <- mapout(testdata)

test_that("mapout() throws error if no rownames exist", {
  d <- testdata
  rownames(d) <- NULL

  expect_error(mapout(d))
})

test_that("predict.mapout() returns character vector", {
  expect_is(predict(obj), "character")
})

obj <- mapout(testdata)
test_that("mapout() returns object of class mapout", {
  expect_is(obj, "mapout")
})

test_that("mapout object has median individual", {
  medians <- plyr::aaply(testdata, 2, median)
  expect_equal(obj$median_array, medians)
})

test_that("mapout object has M and A statistics", {
  expect_is(obj$M, "matrix")
  expect_is(obj$A, "matrix")
})

test_that("dimensions of objects fit", {
  expect_equal(ncol(obj$M), length(obj$median_array))
  expect_equal(ncol(obj$A), length(obj$median_array))

  expect_equal(dim(obj$M), dim(testdata))
  expect_equal(dim(obj$A), dim(testdata))
})

test_that("object has mutual information statistic", {
  expect_is(obj$information, "numeric")
})
