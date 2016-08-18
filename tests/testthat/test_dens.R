testdata <- matrix(c(rnorm(100*100, sd=2), rnorm(100*100, sd=5),
                  rnorm(100*100, sd=3), rnorm(100*100, sd=7)), nrow = 100)
rownames(testdata) <- as.character(1:nrow(testdata))
obj <- dens(testdata)

test_that("dens() returns object of class dens", {
  expect_is(obj, "dens")
})

test_that("plot takes optional highlight= arg of character vector of rownames", {
  pdf(file=NULL) # discard plot
  expect_error(plot(obj, highlight=c(1,2,3)))
  dev.off()
})

test_that("dens() throws error if no rownames exist", {
  d <- testdata
  rownames(d) <- NULL

  expect_error(dens(d))
})
