test_that("mutation adds multivariate normal steps with the given covariance", {
  set.seed(5)
  M <- matrix(c(1, 0.5, 0, 0.5, 1, 0, 0, 0, 0.2), 3)
  model <- evomod_model(additive(1, 3), mutation(rate = 1, cov = M))
  # one locus, rate 1: two mutations (alleles drawn with replacement)
  steps <- replicate(20000, {
    y <- matrix(0, 3, 2)
    rowSums(evomod2:::mutate_cpp(y, model$mutation))
  })
  expect_equal(cov(t(steps)), 2 * M, tolerance = 0.05)
})

test_that("mutational covariance is validated", {
  expect_error(evomod_model(additive(2, 2), mutation(0.1, matrix(c(1, 2, 2, 1), 2))),
               "positive definite")
  expect_error(evomod_model(additive(2, 2), mutation(0.1, diag(3))), "2 x 2")
  expect_equal(evomod_model(additive(2, 2), mutation(0.1, 0.5))$mutation$cov, diag(0.25, 2))
})

test_that("Gaussian log fitness", {
  z <- matrix(c(1, 0, 0, 2, 1, 1), 2)
  omega <- matrix(c(2, 0.5, 0.5, 1), 2)
  theta <- c(0.5, 0.5)
  expected <- apply(z, 2, function(zi) -0.5 * drop(t(zi - theta) %*% solve(omega, zi - theta)))
  expect_equal(drop(evomod2:::log_fitness_cpp(z, theta, omega)), expected)
})

test_that("avg_corr", {
  m <- matrix(c(1, .8, .1, .1, .8, 1, .1, .1, .1, .1, 1, .6, .1, .1, .6, 1), 4)
  a <- avg_corr(m)
  expect_equal(a$within, c(0.8, 0.6))
  expect_equal(a$between, 0.1)
  expect_equal(a$avg_ratio, 7)
})
