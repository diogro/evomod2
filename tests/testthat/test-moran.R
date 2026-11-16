test_that("the Fenwick tree samples in proportion to the weights", {
  set.seed(1)
  w <- c(1, 0, 3, 0.5, 2, 0, 1.5)
  draws <- evomod2:::fenwick_sample_cpp(w, 70000, -1, 0) + 1
  expect_equal(as.numeric(tabulate(draws, 7) / 70000), w / sum(w), tolerance = 0.03)
  expect_false(any(draws %in% c(2, 6)))
  draws <- evomod2:::fenwick_sample_cpp(w, 70000, 2, 0) + 1   # weight 3 set to 0
  w2 <- w; w2[3] <- 0
  expect_equal(as.numeric(tabulate(draws, 7) / 70000), w2 / sum(w2), tolerance = 0.03)
  expect_equal(as.vector(evomod2:::fenwick_sample_cpp(c(0, 0, 5), 10, -1, 0)), rep(2, 10))
})

moran_model <- function(rate = 0, loci = 50, traits = 1, ...) {
  evomod_model(additive(loci, traits), mutation(rate, 0.1), init_sd = 1,
               life_cycle = moran(), ...)
}

test_that("generations overlap: about exp(-1) of the individuals survive a generation", {
  pop <- evolve(moran_model(loci = 5), drift(1), n = 2000, seed = 2)
  before <- apply(genomes(pop), 3, paste, collapse = ",")
  pop <- evolve(pop, drift(1))
  after <- apply(genomes(pop), 3, paste, collapse = ",")
  expect_equal(mean(before %in% after), exp(-1), tolerance = 0.06)
  wf <- evomod_model(additive(5, 1), mutation(0, 0.1), init_sd = 1)
  p <- evolve(wf, drift(1), n = 500, seed = 3)
  b <- apply(genomes(p), 3, paste, collapse = ",")
  expect_equal(mean(b %in% apply(genomes(evolve(p, drift(1))), 3, paste, collapse = ",")), 0)
})

test_that("neutral drift is about twice as fast as in Wright-Fisher", {
  # genic variance of a neutral trait decays by (1 - 1/(2 Ne)) per generation;
  # Ne = n (Wright-Fisher) and about n / 2 (Moran)
  n <- 100; t <- 40
  decay <- function(life_cycle) {
    m <- evomod_model(additive(200, 1), mutation(0, 0.1), init_sd = 1, life_cycle = life_cycle)
    ratios <- replicate(6, {
      pop <- evolve(m, drift(t), n = n, record_every = t)
      g <- g_matrices(pop, "genic")
      g[[2]][1, 1] / g[[1]][1, 1]
    })
    mean(ratios)
  }
  set.seed(4)
  expect_equal(decay(wright_fisher()), (1 - 1 / (2 * n))^t, tolerance = 0.05)
  expect_equal(decay(moran()), (1 - 1 / n)^t, tolerance = 0.07)
})

test_that("Moran populations respond to selection and track the optimum", {
  m <- evomod_model(additive(20, 2), mutation(0.01, 0.1), env_sd = 0.5, init_sd = 0.2,
                    life_cycle = moran())
  pop <- evolve(m, list(stabilizing(20, omega = 5),
                        directional(100, delta_theta = c(0.02, 0), omega = 5)),
                n = 200, record_every = 10, seed = 5)
  last <- trait_means(pop)
  last <- last[last$generation == max(last$generation), ]
  expect_gt(last$mean_x[1] - last$mean_x[2], 1)
  expect_true(all(mean_fitness(pop)$mean_fitness > 0))
  expect_equal(pop$x, genotypic_values(pop))
})

test_that("Moran runs are reproducible and continuing equals one run", {
  m <- evomod_model(binary_b(10, 3), binary_b_mutation(0.01, 0.2, 0.005), env_sd = 0.5,
                    init_sd = 0.2, life_cycle = moran())
  regs <- list(stabilizing(10, omega = 5), fluctuating(20, 5, 0.2, omega = 5))
  one <- evolve(m, regs, n = 60, record_every = 5, seed = 6)
  two <- evolve(evolve(m, regs[1], n = 60, record_every = 5, seed = 6), regs[2])
  expect_identical(one$alleles, two$alleles)
  expect_equal(one$records, two$records)
  expect_equal(one$generation, 30L)
  expect_error(evomod_model(additive(2, 1), mutation(0, 0.1), life_cycle = "moran"), "life_cycle")
})
