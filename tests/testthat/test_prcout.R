testdata <- cbind(rnorm(100, sd=2), rnorm(100, sd=5))
obj <- prcout(testdata)

sdx <- sd(testdata[,1])
sdy <- sd(testdata[,2])
y <- c(2*sdy, -2*sdy)
x <- c(2*sdx, -2*sdx)
outliers <- as.matrix(expand.grid(x,y))
outlierdata <- rbind(outliers, testdata)

test_that("method returns a prcout-object", {
  expect_is(obj, "prcout")
})

test_that("return value contains a prcomp object in $prc", {
  expect_is(obj$prc, "prcomp")
})

test_that("mahalanobis distances exist in $mahalanobis", {
  expect_is(obj$mahalanobis, "numeric")

  # check correct dimension
  expect_true(all(dim(obj$mahalanobis)==c(nrow(testdata), nrow(testdata))))
})

test_that("logical $outliers vector in object", {
  expect_is(obj$outliers, "logical")
})

test_that("prcout() detects some outliers", {
  out <- prcout(outlierdata)
  res <- out$outliers[1:nrow(outliers)]

  expect_true(all(res))
})

test_that("prcout() takes a #sdevs argument for outlier threshold", {
  out1 <- prcout(outlierdata)$outliers[1:nrow(outliers)]
  out2 <- prcout(outlierdata, sdev=5)$outliers[1:nrow(outliers)]

  expect_true(all(out1))
  expect_false(any(out2))
})