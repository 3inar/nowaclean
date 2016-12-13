testdata <- abs(matrix(c(rnorm(100, mean=10, sd=1),
                     rnorm(100, mean=10, sd=2),
                     rnorm(100, mean=10, sd=3),
                     rnorm(100, mean=10, sd=4)),
                     nrow = 100))

negs <- cbind(testdata, matrix(c(rnorm(100, mean=10, sd=1),
                 rnorm(100, mean=10, sd=2)),
                 nrow=100))

colnames(testdata) <- c("V1","V2","V3","V4")
colnames(negs) <- c(colnames(testdata), "V5", "V6")

testdata_with_zeros <- cbind(testdata, rep(0, 100))
colnames(testdata_with_zeros) <- c("V1","V2","V3","V4","zeros")

testdata_with_negative <- cbind(testdata, rnorm(100,sd=1))
colnames(testdata_with_negative) <- c("V1","V2","V3","V4","negative")

test_that("corrected() works when all values are positive", {
  plyr::a_ply(suggested_packages, 1, skip_if_not_installed)
  corrected(testdata, negs)
})

test_that("corrected() throws an error when one column is all zeros", {
  plyr::a_ply(suggested_packages, 1, skip_if_not_installed)
  expect_error(corrected(testdata_with_zeros, negs))
})

test_that("corrected() throws an error when one column has negative values",{
  plyr::a_ply(suggested_packages, 1, skip_if_not_installed)
  expect_error(corrected(testdata_with_negative, negs))
})
