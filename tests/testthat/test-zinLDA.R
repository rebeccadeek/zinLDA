test_that("Alpha is scalar and numeric", {
  set.seed(1)
  d = 30; v=50; k = 3
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10)

  expect_error(zinLDA(sim$sampleTaxaMatrix, K=k, alpha=rep(.1, 3), pi=0.4, a=.5, b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  expect_error(zinLDA(sim$sampleTaxaMatrix, alpha=".1", pi=0.4, a=.5, b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
})

test_that("Pi is scalar and numeric", {
  set.seed(1)
  d = 30; v=50; k = 3
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10)

  expect_error(zinLDA(sim$sampleTaxaMatrix, K=k, alpha=0.1, pi=rep(.4, 3), a=.5, b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  expect_error(zinLDA(sim$sampleTaxaMatrix, K=k, alpha=0.1, pi=".4", a=.5, b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
})

test_that("a is scalar and numeric", {
  set.seed(1)
  d = 30; v=50; k = 3
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10)

  expect_error(zinLDA(sim$sampleTaxaMatrix, K=k, alpha=0.1, pi=0.4, a=rep(.5, 3), b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  expect_error(zinLDA(sim$sampleTaxaMatrix, K=k, alpha=0.1, pi=0.4, a=".5", b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
})

test_that("b is scalar and numeric", {
  set.seed(1)
  d = 30; v=50; k = 3
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10)

  expect_error(zinLDA(sim$sampleTaxaMatrix, K=k, alpha=0.1, pi=0.4, a=.5, b=rep(10, 3)),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  expect_error(zinLDA(sim$sampleTaxaMatrix, K=k, alpha=0.1, pi=0.4, a=.5, b="10"),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
})

test_that("K is a scalar integer", {
  set.seed(1)
  d = 30; v=50
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=3, Alpha=0.1, Pi=0.4, a=.5, b=10)

  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K="3", Alpha=0.1, Pi=0.4, a=.5, b=10))
  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=3.5, Alpha=0.1, Pi=0.4, a=.5, b=10),
               regexp = "ERROR: K must be a length 1 integer.")
})

test_that("sampleTaxaMatrix is numeric matrix", {
  set.seed(1)
  d = 30; v=50; k = 3
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10)

  expect_error(zinLDA(as.character(sim$sampleTaxaMatrix), K=k, alpha=0.1, pi=0.4, a=.5, b=10),
               regexp = 'ERROR: sampleTaxaMatrix must be an numeric matrix.')
  expect_error(zinLDA(as.data.frame(sim$sampleTaxaMatrix), K=k, alpha=0.1, pi=0.4, a=.5, b=10),
               regexp = 'ERROR: sampleTaxaMatrix must be an numeric matrix.')
})

test_that("zinLDA outputs the correct dimensions", {
  set.seed(1)
  d = 20; v=50; k = 3; iter=500
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10)
  modelFit = zinLDA(sim$sampleTaxaMatrix, K=k, alpha=0.1, pi=0.4, a=.5, b=10, runs=iter)

  expect_length(modelFit, 3)
  expect_length(modelFit$delta, iter)
  expect_length(modelFit$beta, iter)
  expect_length(modelFit$theta, iter)

  expect_s3_class(modelFit, "zinLDA_gibbs")
})
