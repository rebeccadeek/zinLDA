test_that("obj is of class 'zinLDA_gibbs'", {
  obj.ls = list()

  expect_error(posterior(obj.ls), regexp = "ERROR: obj must be a list of class 'zinLDA_gibbs'.")
})

test_that("burnin is specified correctly", {
  set.seed(1)
  d = 20; v=50; k = 3; iter=500
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10)
  nMiss =length(which(colSums(sim$sampleTaxaMatrix)==0))
  modelFit = zinLDA(sim$sampleTaxaMatrix, K=k, alpha=0.1, pi=0.4, a=.5, b=10, runs=iter)

  expect_error(posterior(modelFit, burnin = 1.5), regexp = "ERROR: burnin must be an interger.")
  expect_warning(posterior(modelFit, burnin = 0), regexp = "Warning: burnin==0. It is recommend to specify 0 < burnin.")
  expect_error(posterior(modelFit, burnin = 501), regexp = "ERROR: burnin must be smaller than runs.")
})


test_that("zinLDA outputs the correct dimensions", {
  set.seed(1)
  d = 20; v=50; k = 3; iter=500
  N.d = rdu(d, 100, 200)
  sim = simulateZINLDA(D=d, V=v, N=N.d, K=k, Alpha=0.1, Pi=0.4, a=.5, b=10)
  nMiss =length(which(colSums(sim$sampleTaxaMatrix)==0))
  modelFit = zinLDA(sim$sampleTaxaMatrix, K=k, alpha=0.1, pi=0.4, a=.5, b=10, runs=iter)
  post = posterior(modelFit)

  expect_equal(dim(post$theta), c(d,k))
  expect_equal(dim(post$beta), c(k,v-nMiss))
  expect_equal(dim(post$delta), c(k,v-nMiss))
})

