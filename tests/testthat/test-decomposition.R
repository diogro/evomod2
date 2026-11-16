test_that("G = genic + ld, and the genic part is twice the per-locus allele covariance", {
  m <- evomod_model(additive(5, 2), mutation(0.05, 0.2), init_sd = 0.5,
                    recombination = genetic_map(1, 0.1))
  pop <- evolve(m, list(stabilizing(30, omega = 5),
                        directional(30, delta_theta = c(0.05, 0), omega = 5)),
                n = 100, record_every = 10, seed = 1)
  dec <- g_decomposition(pop)
  expect_equal(dec$G, dec$genic + dec$ld)
  expect_equal(unname(dec$G), unname(cov(t(pop$x))))
  # explicit: pool both alleles of every individual at each locus
  g <- genomes(pop)
  n <- dim(g)[3]
  genic <- Reduce(`+`, lapply(1:5, function(l) {
    a <- cbind(g[, 2 * l - 1, ], g[, 2 * l, ])
    2 * cov(t(a)) * (2 * n - 1) / (2 * n)
  })) * n / (n - 1)
  expect_equal(unname(dec$genic), genic)
  # the last record is the current population
  expect_equal(g_matrices(pop, "genic")[["60"]], dec$genic, ignore_attr = TRUE)
  expect_equal(g_matrices(pop, "ld")[["60"]], dec$ld, ignore_attr = TRUE)
  expect_length(g_matrices(pop, "genic"), 7)
})

test_that("binary B decomposition uses allele effects B y", {
  m <- evomod_model(binary_b(5, 2), binary_b_mutation(0.05, 0.2, 0.01), init_sd = 0.5)
  pop <- evolve(m, stabilizing(10, omega = 5), n = 80, seed = 2)
  dec <- g_decomposition(pop)
  g <- genomes(pop)
  eff <- function(i, a) g[-1, a, i] * g[1, a, i]
  genic <- Reduce(`+`, lapply(1:5, function(l) {
    a <- cbind(sapply(1:80, eff, a = 2 * l - 1), sapply(1:80, eff, a = 2 * l))
    2 * cov(t(a)) * (160 - 1) / 160
  })) * 80 / 79
  expect_equal(unname(dec$genic), genic)
  expect_equal(dec$G, dec$genic + dec$ld)
})

test_that("selection builds negative LD under stabilizing selection with linkage", {
  skip_on_cran()
  m <- evomod_model(additive(20, 1), mutation(0.01, 0.2), env_sd = 0.5, init_sd = 0.1,
                    recombination = linkage(rep(0.01, 19)))
  pop <- evolve(m, stabilizing(500, omega = 1), n = 500, record_every = 50, seed = 3)
  ld <- vapply(g_matrices(pop, "ld"), function(x) x[1, 1], 1)
  expect_lt(mean(tail(ld, 5)), 0)   # Bulmer effect: stabilizing selection reduces variance via LD
})

test_that("decomposition is not defined for multilinear epistasis", {
  m <- evomod_model(multilinear(3, 2, random_epsilon(2, 0.1)), mutation(0.05, 0.2), init_sd = 0.2)
  pop <- evolve(m, drift(2), n = 20, seed = 4)
  expect_error(g_decomposition(pop), "sums of allele effects")
  expect_error(g_matrices(pop, "genic"), "sums of allele effects")
  expect_length(g_matrices(pop), 3)
})
