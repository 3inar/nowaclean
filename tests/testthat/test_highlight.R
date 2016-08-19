testdata <- matrix(c(rnorm(100*100, sd=2), rnorm(100*100, sd=5),
                  rnorm(100*100, sd=3), rnorm(100*100, sd=7)), nrow = 100)
rownames(testdata) <- as.character(1:nrow(testdata))

test_that("highlight plot checks input classes", {
  pdf(file=NULL) # drop output

  ma <- mapout(testdata)
  box <- boxout(testdata)
  densty <- dens(testdata)
  pca <- prcout(testdata)

  # should pass:
  highlight("2", pca, box, densty, ma)

  expect_error(mosaic("2", ma, box, densty, ma))
  expect_error(mosaic("2", pca, ma, densty, ma))
  expect_error(mosaic("2", pca, box, ma, ma))
  expect_error(mosaic("2", pca, box, densty, pca))

  expect_error(mosaic(2, pca, box, densty, ma))

  dev.off()
})
