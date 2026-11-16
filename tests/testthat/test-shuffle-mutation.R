test_that("shuffling keeps each locus's allele pool and the allele states", {
  m <- evomod_model(binary_b(6, 3), binary_b_mutation(0.05, 0.2, 0.05), init_sd = 0.3)
  pop <- evolve(m, stabilizing(10, omega = 5), n = 40, seed = 1)
  g <- genomes(pop)
  s <- shuffle_genomes(g)
  expect_equal(dim(s), dim(g))
  key <- function(a, l) sort(apply(matrix(a[, c(2 * l - 1, 2 * l), ], nrow = dim(a)[1]), 2,
                                   paste, collapse = ","))
  for (l in 1:6) expect_identical(key(s, l), key(g, l))
  expect_false(identical(s, g))
})

test_that("shuffling removes linkage disequilibrium", {
  n <- 2000
  g <- array(0, c(1, 4, n))
  g[, c(1, 3), ] <- 1                 # paternal haplotype (1, 1), maternal (0, 0)
  set.seed(2)
  s <- shuffle_genomes(g)
  D <- function(a) { x <- c(a[1, 1, ], a[1, 2, ]); y <- c(a[1, 3, ], a[1, 4, ]); mean(x * y) - mean(x) * mean(y) }
  expect_equal(D(g), 0.25)
  expect_lt(abs(D(s)), 0.02)
})

test_that("for additive architectures the shuffle estimate matches the exact genic part", {
  m <- evomod_model(additive(10, 2), mutation(0.02, 0.3), init_sd = 0.4,
                    recombination = linkage(rep(0.01, 9)))
  pop <- evolve(m, stabilizing(200, omega = 2), n = 200, seed = 3)
  exact <- g_decomposition(pop)
  set.seed(4)
  shuf <- g_decomposition(pop, method = "shuffle", shuffles = 400)
  expect_equal(shuf$G, exact$G)
  expect_equal(shuf$genic, exact$genic, tolerance = 0.05)
  expect_equal(shuf$G, shuf$genic + shuf$ld)
})

test_that("shuffle decomposition works for multilinear epistasis", {
  m <- evomod_model(multilinear(8, 2, random_epsilon(2, 0.2)), mutation(0.02, 0.2), init_sd = 0.3,
                    recombination = linkage(rep(0.01, 7)))
  pop <- evolve(m, stabilizing(100, omega = 3), n = 100, seed = 5)
  expect_error(g_decomposition(pop), "shuffle")
  d <- g_decomposition(pop, method = "shuffle", shuffles = 5)
  expect_equal(d$G, d$genic + d$ld)
  expect_true(all(diag(d$genic) > 0))
})

test_that("single mutations in an additive model have the mutational covariance", {
  M <- matrix(c(0.04, 0.02, 0.02, 0.09), 2)
  m <- evomod_model(additive(10, 2), mutation(rate = 0.001, cov = M), init_sd = 0.1)
  pop <- evolve(m, drift(2), n = 50, seed = 6)
  set.seed(7)
  mc <- mutational_covariance(pop, n = 50000)
  expect_equal(unname(mc$per_mutation), M, tolerance = 0.03)
  expect_equal(unname(mc$mean), c(0, 0), tolerance = 0.005)
  expect_equal(mc$expected_mutations, 0.001 * 20)
  expect_equal(unname(mc$M), 0.02 * unname(mc$per_mutation), tolerance = 0.03)
})

test_that("multilinear mutation effects scale with the genetic background", {
  # every allele of every individual has effect v: S - Y_l = 2 v (L - 1), so a
  # mutation of size d at any locus changes x by d (1 + eps 2 v (L - 1))
  L <- 10; v <- 0.1; eps <- 0.3; s <- 0.05
  m <- evomod_model(multilinear(L, 1, array(eps, c(1, 1, 1))), mutation(rate = 0.01, cov = s))
  pop <- evolve(m, drift(1), n = 20, seed = 8, initial = array(v, c(1, 2 * L, 20)))
  pop$alleles[] <- v                       # undo the generation's mutations
  set.seed(9)
  mc <- mutational_covariance(pop, n = 40000, effects = TRUE)
  expect_equal(unname(mc$per_mutation[1, 1]), s^2 * (1 + eps * 2 * v * (L - 1))^2, tolerance = 0.03)
  add <- evomod_model(additive(L, 1), mutation(rate = 0.01, cov = s))
  pop_add <- pop; pop_add$model <- add
  mc_add <- mutational_covariance(pop_add, n = 40000)
  expect_gt(mc$per_mutation[1, 1], 1.5 * mc_add$per_mutation[1, 1])
})

test_that("binary B mutation effects follow B y", {
  # all B entries 1 and y = 0.2: a y mutation of size d changes every trait by d;
  # a B flip changes one trait by -0.2
  m <- evomod_model(binary_b(5, 2), binary_b_mutation(rate = 0.01, sd = 0.1, rate_b = 0.01,
                                                      b_mutation = "flip"))
  g <- array(1, c(3, 10, 30)); g[1, , ] <- 0.2
  pop <- evolve(m, drift(1), n = 30, seed = 10, initial = g)
  pop$alleles[] <- g
  set.seed(11)
  mc <- mutational_covariance(pop, n = 30000, effects = TRUE)
  dx <- mc$effects
  y_mut <- abs(dx[1, ] - dx[2, ]) < 1e-12 & dx[1, ] != 0
  flips <- colSums(dx != 0) == 1
  expect_equal(mean(y_mut), 1 / 3, tolerance = 0.05)     # y events: 10 sites vs 20 B entries
  expect_equal(mean(flips), 2 / 3, tolerance = 0.05)
  expect_true(all(abs(dx[dx != 0 & rep(flips, each = 2)] + 0.2) < 1e-12))
  expect_equal(sd(dx[1, y_mut]), 0.1, tolerance = 0.05)
  expect_equal(mc$expected_mutations, 0.01 * 10 + 0.01 * 20)
})

test_that("local effects work for single-trait populations", {
  m <- evomod_model(additive(4, 1), mutation(0.01, 0.1), init_sd = 0.1)
  pop <- evolve(m, drift(2), n = 10, seed = 12)
  expect_equal(dim(local_effects(pop)), c(1L, 4L))
})
