naive_multilinear <- function(alleles, epsilon) {
  p <- nrow(alleles)
  Y <- alleles[, c(TRUE, FALSE), drop = FALSE] + alleles[, c(FALSE, TRUE), drop = FALSE]
  x <- rowSums(Y)
  for (i in seq_len(ncol(Y))) for (j in seq_len(ncol(Y))) if (i < j)
    for (k in seq_len(p)) x[k] <- x[k] + drop(t(Y[, i]) %*% epsilon[k, , ] %*% Y[, j])
  x
}

test_that("architectures validate their arguments", {
  expect_error(additive(0, 2), "loci")
  expect_error(multilinear(5, 2, array(1:8, c(2, 2, 2))), "symmetric")
  expect_error(multilinear(5, 2, array(0, c(2, 2, 3))), "2 x 2 x 2")
  e <- random_epsilon(3, sd = 0.5)
  for (k in 1:3) expect_true(isSymmetric(e[k, , ]))
  expect_s3_class(multilinear(5, 3, e), "evomod_architecture")
})

test_that("additive genotypic value is the sum of allele effects", {
  set.seed(1)
  y <- matrix(rnorm(3 * 12), 3)
  expect_equal(drop(genotypic_values(additive(6, 3), y)), rowSums(y))
  expect_equal(genotypic_values(multilinear(6, 3, array(0, c(3, 3, 3))), y),
               genotypic_values(additive(6, 3), y))
})

test_that("multilinear genotypic value equals the explicit pairwise sum", {
  set.seed(2)
  e <- random_epsilon(3, sd = 0.5)
  y <- matrix(rnorm(3 * 12), 3)
  expect_equal(drop(genotypic_values(multilinear(6, 3, e), y)), naive_multilinear(y, e))
  # several genomes at once
  ys <- array(rnorm(3 * 12 * 4), c(3, 12, 4))
  expect_equal(genotypic_values(multilinear(6, 3, e), ys)[, 3], naive_multilinear(ys[, , 3], e))
})

test_that("univariate multilinear model matches the textbook formula", {
  set.seed(3)
  y <- matrix(rnorm(12), 1)
  Y <- y[c(TRUE, FALSE)] + y[c(FALSE, TRUE)]
  pairs <- outer(Y, Y)[upper.tri(diag(6))]
  expect_equal(drop(genotypic_values(multilinear(6, 1, array(0.3, c(1, 1, 1))), y)),
               sum(Y) + 0.3 * sum(pairs))
})

test_that("local effects are the change from removing a locus", {
  set.seed(4)
  e <- random_epsilon(3, sd = 0.5)
  arch <- multilinear(6, 3, e)
  y <- matrix(rnorm(3 * 12), 3)
  y0 <- y
  y0[, 5:6] <- 0
  L <- local_effects(arch, y)
  expect_equal(drop(genotypic_values(arch, y) - genotypic_values(arch, y0)), L[, 3])
  expect_equal(local_effects(additive(6, 3), y), y[, c(TRUE, FALSE)] + y[, c(FALSE, TRUE)])
})
