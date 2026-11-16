small_model <- function(arch = additive(10, 2), rate = 0.01) {
  evomod_model(arch, mutation(rate = rate, cov = 0.1), env_sd = 0.5, init_sd = 0.2)
}

test_that("runs are reproducible with a seed", {
  m <- small_model()
  a <- evolve(m, list(stabilizing(20, omega = 5)), n = 50, seed = 1)
  b <- evolve(m, list(stabilizing(20, omega = 5)), n = 50, seed = 1)
  expect_identical(a$alleles, b$alleles)
  expect_identical(a$records, b$records)
})

test_that("continuing a population equals one run", {
  m <- small_model(multilinear(10, 2, random_epsilon(2, sd = 0.1)))
  regs <- list(stabilizing(15, omega = 5), directional(15, delta_theta = c(0.1, 0), omega = 5))
  one <- evolve(m, regs, n = 50, record_every = 5, seed = 2)
  two <- evolve(m, regs[1], n = 50, record_every = 5, seed = 2)
  two <- evolve(two, regs[2])
  expect_identical(one$alleles, two$alleles)
  expect_equal(one$records, two$records)
  expect_equal(one$generation, 30L)
  expect_equal(two$generation, 30L)
})

test_that("records and the optimum follow the regimes", {
  m <- small_model()
  pop <- evolve(m, list(drift(10), directional(10, delta_theta = c(0.5, -0.5), omega = 5),
                        stabilizing(10, omega = 5)),
                n = 40, record_every = 5, seed = 3)
  tm <- trait_means(pop)
  expect_equal(unique(tm$generation), seq(0, 30, by = 5))
  th <- tm$theta[tm$trait == "trait1"]
  expect_equal(th, c(0, 0, 0, 2.5, 5, 5, 5))
  expect_equal(pop$theta, c(5, -5))
  expect_length(g_matrices(pop), 7)
  expect_equal(dim(g_matrices(pop)[["30"]]), c(2, 2))
  v <- variances(pop)
  expect_equal(v$G, unlist(lapply(g_matrices(pop), diag)), ignore_attr = TRUE)
  expect_true(all(v$h2 > 0 & v$h2 < 1))
  # explicit theta restarts the optimum
  pop <- evolve(pop, stabilizing(5, omega = 5, theta = c(1, 1)))
  expect_equal(pop$theta, c(1, 1))
})

test_that("genotypic values of the population match the architecture", {
  m <- small_model(multilinear(10, 2, random_epsilon(2, sd = 0.1)))
  pop <- evolve(m, stabilizing(5, omega = 5), n = 30, seed = 4)
  expect_equal(pop$x, genotypic_values(pop))
  expect_equal(dim(genomes(pop)), c(2L, 20L, 30L))
})

test_that("without mutation, offspring alleles come from the same locus of the parents", {
  m <- small_model(rate = 0)
  pop <- evolve(m, drift(1), n = 30, seed = 5)
  before <- genomes(pop)
  pop <- evolve(pop, drift(1))
  after <- genomes(pop)
  for (l in 1:10) {
    pool <- apply(before[, c(2 * l - 1, 2 * l), ], c(2, 3), paste, collapse = ",")
    kids <- apply(after[, c(2 * l - 1, 2 * l), ], c(2, 3), paste, collapse = ",")
    expect_true(all(kids %in% pool))
  }
})

test_that("directional selection moves the population towards the optimum", {
  m <- small_model()
  pop <- evolve(m, list(stabilizing(10, omega = 5), directional(100, delta_theta = c(0.02, 0), omega = 5)),
                n = 200, record_every = 10, seed = 6)
  last <- trait_means(pop)
  last <- last[last$generation == max(last$generation), ]
  expect_gt(last$mean_x[1] - last$mean_x[2], 1)
})

test_that("directional epistasis changes additive variance as theory predicts", {
  skip_on_cran()
  ratio <- function(e) {
    E <- array(0, c(2, 2, 2)); E[1, 1, 1] <- e; E[2, 2, 2] <- e
    m <- evomod_model(multilinear(50, 2, E), mutation(rate = 0.005, cov = 0.05),
                      env_sd = 1, init_sd = 0.05)
    pop <- evolve(m, list(stabilizing(1000, omega = 20),
                          directional(1000, delta_theta = c(0.005, 0), omega = 20)),
                  n = 500, record_every = 1000, seed = 7)
    G <- variances(pop)
    G <- G[G$trait == "trait1", "G"]
    G[3] / G[2]
  }
  expect_gt(ratio(0.05), ratio(0))
  expect_lt(ratio(-0.05), ratio(0))
})

test_that("input errors", {
  m <- small_model()
  expect_error(evolve(m, list("drift"), n = 10), "regimes")
  expect_error(evolve(m, directional(5, delta_theta = 1, omega = 5), n = 10), "length 2")
  expect_error(evolve(m, stabilizing(5, omega = diag(3)), n = 10), "2 x 2")
})
