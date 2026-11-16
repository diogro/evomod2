cm_model <- function(rate_w = 0.01, rate_c = 0.01, sd = 0.5, ...) {
  evomod_model(contribution_matrix(traits = 3, genes = 2, ...),
               contribution_mutation(rate_w, rate_c, sd), env_sd = 0.2)
}
# genome (1 x alleles) from weights a and contributions C (each allele half)
cm_genome <- function(a, C) matrix(rep(c(a, as.vector(C)) / 2, each = 2), nrow = 1)

test_that("genotypic values are C a", {
  arch <- contribution_matrix(traits = 3, genes = 2)
  expect_equal(arch$loci, 2 + 3 * 2)
  a <- c(0.5, -2); C <- matrix(c(1, 0.2, -1, 0, 3, 0.5), 3)
  g <- cm_genome(a, C)
  expect_equal(drop(genotypic_values(arch, g)), drop(C %*% a))
  L <- local_effects(arch, g)
  expect_equal(dim(L), c(3L, 8L))
  for (l in 1:8) {
    g0 <- g; g0[, c(2 * l - 1, 2 * l)] <- 0
    expect_equal(drop(genotypic_values(arch, g) - genotypic_values(arch, g0)), L[, l])
  }
  expect_equal(contribution_matrix(3, 2)$init_contributions, diag(1, 3, 2))
  expect_error(contribution_matrix(3, 2, init_contributions = diag(3)), "3 x 2")
  expect_error(evomod_model(contribution_matrix(2), mutation(0.1, 0.1)), "contribution_mutation")
})

test_that("new populations start at the initial weights and contributions", {
  m <- evomod_model(contribution_matrix(4), contribution_mutation(0, 0, 0.5))
  pop <- evolve(m, drift(1), n = 10, seed = 1)
  expect_equal(pop$x, matrix(1, 4, 10))
  m2 <- evomod_model(contribution_matrix(2, 2, init_weights = c(2, 1),
                                         init_contributions = matrix(c(1, 1, 0, 1), 2)),
                     contribution_mutation(0, 0, 0.5))
  expect_equal(evolve(m2, drift(1), n = 3, seed = 2)$x[, 1], c(2, 3))
})

test_that("mutation rates differ between weight and contribution genes", {
  m <- cm_model(rate_w = 0.2, rate_c = 0.02)
  g <- array(0, c(1, 16, 1))
  set.seed(3)
  hits <- replicate(5000, {
    out <- evomod2:::mutate_cpp(matrix(g, 1), m$mutation)
    c(weight = sum(out[1:4] != 0), contribution = sum(out[5:16] != 0))
  })
  # mutations pick alleles with replacement: expected number of distinct
  # alleles hit by Binomial(n, rate) mutations
  distinct <- function(n, rate) sum(dbinom(0:n, n, rate) * n * (1 - (1 - 1 / n)^(0:n)))
  expect_equal(mean(hits["weight", ]), distinct(4, 0.2), tolerance = 0.05)
  expect_equal(mean(hits["contribution", ]), distinct(12, 0.02), tolerance = 0.1)
  expect_equal(evomod2:::expected_mutations_cpp(matrix(g, 1), m$mutation), 4 * 0.2 + 12 * 0.02)
})

test_that("a contribution mutation changes one trait, a weight mutation a column of C", {
  m <- cm_model(rate_w = 0.01, rate_c = 0.01, sd = 0.3)
  C <- matrix(c(1, 0.5, 0, 0, 1, 2), 3)
  pop <- evolve(m, drift(1), n = 10, seed = 4, initial = array(cm_genome(c(1, 1), C), c(1, 16, 10)))
  pop$alleles[] <- cm_genome(c(1, 1), C)
  set.seed(5)
  mc <- mutational_covariance(pop, n = 20000, effects = TRUE)
  dx <- mc$effects
  changed <- colSums(abs(dx) > 1e-12)
  weight <- apply(dx, 2, function(d) any(abs(d) > 1e-12) &&
                    (all(abs(d - C[, 1] * d[1]) < 1e-9) || all(abs(d - C[, 2] * d[3] / 2) < 1e-9)))
  expect_equal(mean(changed == 1), 12 / 16, tolerance = 0.05)   # 12 contribution alleles of 16
  expect_true(all(changed <= 3))
  expect_equal(mc$expected_mutations, 16 * 0.01)
})

test_that("the shuffle decomposition and selection work", {
  m <- evomod_model(contribution_matrix(4), contribution_mutation(0.01, 0.005, 0.3),
                    env_sd = 0.5, init_sd = 0.05)
  pop <- evolve(m, list(drift(10), fluctuating(300, 10, 0.2, omega = 1)), n = 200,
                record_every = 10, seed = 6)
  expect_error(g_decomposition(pop), "shuffle")
  d <- g_decomposition(pop, method = "shuffle", shuffles = 5)
  expect_equal(d$G, d$genic + d$ld)
  w <- mean_fitness(pop)
  expect_true(all(is.na(w$mean_fitness[w$generation <= 10])))
  expect_true(all(w$mean_fitness[w$generation > 10] > 0 & w$mean_fitness[w$generation > 10] <= 1))
  # moving away from the start (z = 1) towards the optimum (0 to +-1)
  expect_lt(mean(abs(pop$x[, 1] - pop$theta)), 1)
})
