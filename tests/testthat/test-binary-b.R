bb_model <- function(rate = 0.01, rate_b = 0.001, b_mutation = "random", ...) {
  evomod_model(binary_b(10, 4, ...), binary_b_mutation(rate, 0.2, rate_b, b_mutation),
               env_sd = 0.5, init_sd = 0.2)
}

test_that("binary B genotypic value is B y", {
  set.seed(1)
  g <- rbind(rnorm(20), matrix(rbinom(80, 1, 0.5), 4))
  expect_equal(drop(genotypic_values(binary_b(10, 4), g)), drop(g[-1, ] %*% g[1, ]))
  L <- local_effects(binary_b(10, 4), g)
  expect_equal(dim(L), c(4L, 10L))
  expect_equal(L[, 3], g[-1, 5] * g[1, 5] + g[-1, 6] * g[1, 6])
  expect_equal(rowSums(L), drop(g[-1, ] %*% g[1, ]))
  expect_error(genotypic_values(binary_b(10, 4), g[-1, ]), "5 x 20")
})

test_that("architecture and mutation model must match", {
  expect_error(evomod_model(binary_b(5, 2), mutation(0.1, 0.1)), "binary_b_mutation")
  expect_error(evomod_model(additive(5, 2), binary_b_mutation(0.1, 0.1, 0.1)), "mutation\\(\\)")
  expect_error(binary_b(5, 2, init_b = 2), "probability")
})

test_that("new binary B populations follow init_b", {
  pop <- evolve(bb_model(init_b = 0.25), drift(1), n = 200, seed = 2)
  B <- genomes(pop)[-1, , ]
  expect_true(all(B %in% c(0, 1)))
  expect_equal(mean(B), 0.25, tolerance = 0.05)
})

test_that("binary B mutation changes y by Normal(0, sd) and B entries at rate_b", {
  set.seed(3)
  g <- rbind(rep(0, 2000), matrix(0, 4, 2000))
  m <- binary_b_mutation(rate = 0.1, sd = 0.5, rate_b = 0.1, b_mutation = "flip")
  out <- evomod2:::mutate_cpp(g, m)
  dy <- out[1, ] - g[1, ]
  expect_equal(mean(dy != 0), 0.1, tolerance = 0.15)
  expect_equal(sd(dy[dy != 0]), 0.5, tolerance = 0.1)
  expect_equal(mean(out[-1, ]), 0.1, tolerance = 0.1)   # every mutated entry flips
  random <- evomod2:::mutate_cpp(g, binary_b_mutation(0.1, 0.5, 0.1, "random"))
  expect_equal(mean(random[-1, ]), 0.05, tolerance = 0.15)  # half of them change
  expect_true(all(random[-1, ] %in% c(0, 1)))
})

test_that("binary B loci are inherited intact and B stays binary", {
  m <- bb_model(rate = 0, rate_b = 0)
  pop <- evolve(m, drift(1), n = 30, seed = 4)
  before <- genomes(pop)
  pop <- evolve(pop, drift(1))
  after <- genomes(pop)
  for (l in 1:10) {
    pool <- apply(before[, c(2 * l - 1, 2 * l), ], c(2, 3), paste, collapse = ",")
    kids <- apply(after[, c(2 * l - 1, 2 * l), ], c(2, 3), paste, collapse = ",")
    expect_true(all(kids %in% pool))
  }
  pop <- evolve(bb_model(rate_b = 0.05), stabilizing(20, omega = 5), n = 30, seed = 6)
  expect_true(all(genomes(pop)[-1, , ] %in% c(0, 1)))
  expect_equal(pop$x, genotypic_values(pop))
})

test_that("binary B populations respond to directional selection", {
  pop <- evolve(bb_model(), list(stabilizing(10, omega = 5),
                                 directional(100, delta_theta = c(0.02, 0.02, 0, 0), omega = 5)),
                n = 200, record_every = 10, seed = 5)
  last <- trait_means(pop)
  last <- last[last$generation == max(last$generation), ]
  expect_gt(mean(last$mean_x[1:2]) - mean(last$mean_x[3:4]), 1)
  expect_equal(nrow(variances(pop)), 4 * 12)
})

test_that("C-style initial effects are sums of random mutations", {
  m <- evomod_model(binary_b(50, 2, init_y = "mutations"), binary_b_mutation(0.01, 0.2, 0.001),
                    init_sd = 0.2)
  # about exp(-1) of the alleles are never hit
  y <- genomes(evolve(m, drift(1), n = 200, seed = 7))[1, , ]
  expect_true(mean(y == 0) > 0.25)
  expect_equal(sd(y), 0.2, tolerance = 0.1)
})
