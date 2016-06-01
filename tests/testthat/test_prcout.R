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

test_that("predict() detects some outliers", {
  out <- prcout(outlierdata)
  res <- predict(out)[1:nrow(outliers)]

  expect_true(all(res))
})

test_that("predict() takes optional #sdev argument for outlier threshold", {
  object <- prcout(outlierdata)
  out1 <- predict(object)[1:nrow(outliers)]
  out2 <- predict(object, sdev=5)[1:nrow(outliers)]

  expect_true(all(out1))
  expect_false(all(out2))
})

test_that("prcout() takes optional 'prob' argument to define the mass of data", {
  out1 <- predict(prcout(outlierdata))
  out2 <- predict(prcout(outlierdata, prob=0.5))

  expect_more_than(sum(out2), sum(out1))
})
