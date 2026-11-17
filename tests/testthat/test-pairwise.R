naive_pairwise <- function(y, eps) {
  p <- nrow(y); L <- ncol(y) / 2
  Y <- y[, c(TRUE, FALSE), drop = FALSE] + y[, c(FALSE, TRUE), drop = FALSE]
  x <- rowSums(Y)
  for (i in 1:(L - 1)) for (j in (i + 1):L) for (k in 1:p)
    x[k] <- x[k] + drop(t(Y[, i]) %*% eps[k, , , i, j] %*% Y[, j])
  x
}

test_that("pairwise genotypic values equal the explicit sum over pairs of loci", {
  set.seed(1)
  eps <- pairwise_epsilon(6, 3, sd = 0.5)
  arch <- multilinear(6, 3, eps)
  expect_equal(arch$type, "multilinear_pairwise")
  y <- matrix(rnorm(3 * 12), 3)
  expect_equal(drop(genotypic_values(arch, y)), naive_pairwise(y, eps))
  # asymmetric coefficients are allowed: (i, j) and (j, i) orders differ
  eps[1, 1, 2, 1, 2] <- 3
  expect_equal(drop(genotypic_values(multilinear(6, 3, eps), y)), naive_pairwise(y, eps))
  expect_error(multilinear(6, 3, array(0, c(3, 3, 3, 5, 5))), "3 x 3 x 3 x 6 x 6")
})

test_that("the same tensor for every pair gives the shared multilinear model", {
  set.seed(2)
  e <- random_epsilon(2, sd = 0.4)
  eps <- array(0, c(2, 2, 2, 5, 5))
  for (i in 1:4) for (j in (i + 1):5) eps[, , , i, j] <- e
  y <- matrix(rnorm(20), 2)
  expect_equal(genotypic_values(multilinear(5, 2, eps), y), genotypic_values(multilinear(5, 2, e), y))
  expect_equal(local_effects(multilinear(5, 2, eps), y), local_effects(multilinear(5, 2, e), y))
})

test_that("pairwise local effects are the change from removing a locus", {
  set.seed(3)
  eps <- pairwise_epsilon(5, 2, sd = 1)
  arch <- multilinear(5, 2, eps)
  y <- matrix(rnorm(20), 2)
  L <- local_effects(arch, y)
  for (l in 1:5) {
    y0 <- y; y0[, c(2 * l - 1, 2 * l)] <- 0
    expect_equal(drop(genotypic_values(arch, y) - genotypic_values(arch, y0)), L[, l])
  }
})

test_that("pairwise_epsilon options", {
  set.seed(4)
  e <- pairwise_epsilon(4, 2, sd = 1, within_trait = TRUE, traits_with_epistasis = 1)
  expect_true(all(e[1, 1, 1, 1, 2:4] != 0))
  e[1, 1, 1, , ] <- 0
  expect_true(all(e == 0))
  e <- pairwise_epsilon(4, 2, sd = 2)
  expect_true(all(e[, , , 2, 1] == 0) && all(e[, , , 3, 3] == 0))
  expect_equal(sd(e[e != 0]), 2, tolerance = 0.3)
})

test_that("additive covariance: all of G when additive, less with epistasis", {
  m <- evomod_model(additive(10, 2), mutation(0.01, 0.2), env_sd = 1, init_sd = 0.3)
  pop <- evolve(m, stabilizing(50, omega = 20), n = 200, seed = 5)
  a <- additive_covariance(pop)
  expect_equal(a$additive, a$G, tolerance = 1e-8)
  expect_equal(unname(a$nonadditive), matrix(0, 2, 2), tolerance = 1e-8)
  set.seed(6)
  m <- evomod_model(multilinear(10, 2, pairwise_epsilon(10, 2, sd = 2)), mutation(0.01, 0.2),
                    env_sd = 1, init_sd = 0.3)
  pop <- evolve(m, stabilizing(50, omega = 20), n = 300, seed = 7)
  a <- additive_covariance(pop)
  expect_true(all(diag(a$additive) > 0) && all(diag(a$additive) < diag(a$G)))
  expect_equal(a$G, a$additive + a$nonadditive)
  # binary B: additive, so all of G
  mb <- evomod_model(binary_b(8, 2), binary_b_mutation(0.02, 0.3, 0.01), init_sd = 0.3)
  a <- additive_covariance(evolve(mb, stabilizing(30, omega = 20), n = 200, seed = 8))
  expect_equal(a$additive, a$G, tolerance = 1e-8)
  expect_error(additive_covariance(evolve(m, drift(1), n = 15, seed = 9)), "more individuals")
})

test_that("pairwise models evolve and have mutational covariances", {
  set.seed(10)
  m <- evomod_model(multilinear(8, 2, pairwise_epsilon(8, 2, sd = 1)), mutation(0.01, 0.2),
                    env_sd = 1)
  pop <- evolve(m, stabilizing(100, omega = 20), n = 100, seed = 11)
  expect_equal(pop$x, genotypic_values(pop))
  mc <- mutational_covariance(pop, n = 2000)
  expect_equal(dim(mc$per_mutation), c(2L, 2L))
  expect_error(g_decomposition(pop), "shuffle")
})
