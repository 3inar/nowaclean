testdata <- cbind(rnorm(100, sd=2), rnorm(100, sd=5))
rownames(testdata) <- as.character(1:nrow(testdata))
obj <- prcout(testdata)

sdx <- sd(testdata[,1])
sdy <- sd(testdata[,2])
y <- c(3*sdy, -3*sdy)
x <- c(3*sdx, -3*sdx)
outliers <- as.matrix(expand.grid(x,y))
outlierdata <- rbind(outliers, testdata)
rownames(outlierdata) <- as.character(1:nrow(outlierdata))

test_that("method returns a prcout-object", {
  expect_is(obj, "prcout")
})

test_that("prcout() throws error if no rownames are provided", {
  expect_error(prcout(cbind(rnorm(100, sd=2), rnorm(100, sd=5))))
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

  expect_true(length(res) > 0)
})

test_that("predict() takes optional #sdev argument for outlier threshold", {
  object <- prcout(outlierdata)
  predict(object)[1:nrow(outliers)]
  out1 <- predict(object, sdev=1)
  out2 <- predict(object, sdev=6)

  expect_true(length(out1) >= nrow(outliers))
})

test_that("prcout() takes optional 'prob' argument to define the mass of data", {
  out1 <- predict(prcout(outlierdata))
  out2 <- predict(prcout(outlierdata, prob=0.5))

  expect_true(length(out2) > length(out1))
})
