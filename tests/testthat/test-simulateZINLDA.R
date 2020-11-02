test_that("Alpha is scalar and numeric", {
  set.seed(1)
  d = 30; v=50; k=3
  N.d = rdu(d, 100, 200)

  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=rep(.1, 3), Pi=0.4, a=.5, b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=".1", Pi=0.4, a=.5, b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
})

test_that("Pi is scalar and numeric", {
  set.seed(1)
  d = 30; v=50; k=3
  N.d = rdu(d, 100, 200)

  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=rep(.4, 3), a=.5, b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=".4", a=.5, b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
})

test_that("a is scalar and numeric", {
  set.seed(1)
  d = 30; v=50; k=3
  N.d = rdu(d, 100, 200)

  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=rep(.5, 3), b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=".5", b=10),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
})

test_that("b is scalar and numeric", {
  set.seed(1)
  d = 30; v=50; k=3
  N.d = rdu(d, 100, 200)

  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=rep(10, 3)),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b="10"),
               regexp = "ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
})

test_that("N is a length D numeric vector", {
  set.seed(1)
  d = 30; v=50; k=3
  N.d = rdu(d+1, 100, 200)

  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10),
               regexp = 'ERROR: N must be a length D numeric vector.')
  expect_error(simulateZINLDA(D=d, V=v, N=rep("d", d), K=k, Alpha=0.1, Pi=0.4, a=.5, b=10),
               regexp = 'ERROR: N must be a length D numeric vector.')
})

test_that("D is a scalar integer", {
  set.seed(1)
  v=50; k=3
  N.d = rdu(30, 100, 200)

  expect_error(simulateZINLDA(D="30", V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10))
  expect_error(simulateZINLDA(D=30.5, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10),
               regexp = "ERROR: D must be a length 1 integer.")
})

test_that("V is a scalar integer", {
  set.seed(1)
  d=30; k=3
  N.d = rdu(30, 100, 200)

  expect_error(simulateZINLDA(D=d, V="50", N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10))
  expect_error(simulateZINLDA(D=d, V=50.5, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10),
               regexp = "ERROR: V must be a length 1 integer.")
})

test_that("K is a scalar integer", {
  set.seed(1)
  d=30; v=50
  N.d = rdu(30, 100, 200)

  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K="3", Alpha=0.1, Pi=0.4, a=.5, b=10))
  expect_error(simulateZINLDA(D=d, V=v, N=N.d, K=3.5, Alpha=0.1, Pi=0.4, a=.5, b=10),
               regexp = "ERROR: K must be a length 1 integer.")
})

test_that("simulateZINLDA outputs the correct dimensions", {
  set.seed(1)
  d = 30; v=50; k=3
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10)

  expect_length(sim, 6)
  expect_length(sim$cohort, d)
  expect_length(sim$z, d)

  expect_equal(dim(sim$sampleTaxaMatrix), c(d,v))
  expect_equal(dim(sim$theta), c(d,k))
  expect_equal(dim(sim$beta), c(k,v))
  expect_equal(dim(sim$delta), c(k,v))
})
