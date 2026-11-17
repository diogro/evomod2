test_that("viability selection: neutral drift with an effective size of about n", {
  n <- 100; t <- 40
  m <- evomod_model(additive(200, 1), mutation(0, 0.1), init_sd = 1,
                    life_cycle = viability_selection())
  set.seed(1)
  ratio <- mean(replicate(8, {
    pop <- evolve(m, drift(t), n = n, record_every = t)
    g <- g_matrices(pop, "genic")
    g[[2]][1, 1] / g[[1]][1, 1]
  }))
  expect_equal(ratio, (1 - 1 / (2 * n))^t, tolerance = 0.05)
})

test_that("offspring come from a female (first half) and a male (second half)", {
  # every individual carries a unique paternal value at locus 1 and maternal
  # value at locus 2, so a child's parents can be identified
  n <- 40
  g <- array(0, c(1, 4, n))
  g[1, 1:2, ] <- rep(seq_len(n), each = 2)       # locus 1: individual id
  m <- evomod_model(additive(2, 1), mutation(0, 0.1), life_cycle = viability_selection())
  pop <- evolve(m, drift(1), n = n, initial = g, seed = 2)
  kids <- genomes(pop)
  fathers <- kids[1, 1, ]                         # paternal allele at locus 1
  mothers <- kids[1, 2, ]                         # maternal allele at locus 1
  expect_true(all(fathers > n / 2))
  expect_true(all(mothers <= n / 2))
})

test_that("viability selection responds to selection and keeps n under strong selection", {
  m <- evomod_model(additive(20, 2), mutation(0.01, 0.1), env_sd = 0.5, init_sd = 0.2,
                    life_cycle = viability_selection())
  pop <- evolve(m, list(stabilizing(20, omega = 5),
                        directional(100, delta_theta = c(0.02, 0), omega = 5)),
                n = 200, record_every = 10, seed = 3)
  last <- trait_means(pop)
  last <- last[last$generation == max(last$generation), ]
  expect_gt(last$mean_x[1] - last$mean_x[2], 1)
  # strong selection: fewer survivors than n, population size kept by resampling
  strong <- evolve(m, stabilizing(3, omega = 0.05, theta = c(0.3, 0.3)), n = 50, seed = 4)
  expect_equal(dim(genomes(strong))[3], 50L)
  # no survivors at all: extinction is an error
  expect_error(evolve(m, stabilizing(1, omega = 0.0001, theta = c(50, 50)), n = 50, seed = 4),
               "no zygote survived")
})

test_that("viability selection is reproducible, continues exactly and checks its input", {
  m <- evomod_model(multilinear(8, 2, pairwise_epsilon(8, 2, sd = 0.5)), mutation(0.01, 0.1),
                    env_sd = 1, life_cycle = viability_selection(mates = 3, offspring_per_mate = 1))
  regs <- list(stabilizing(10, omega = 10), stabilizing(10, omega = 10))
  one <- evolve(m, regs, n = 60, record_every = 5, seed = 5)
  two <- evolve(evolve(m, regs[1], n = 60, record_every = 5, seed = 5), regs[2])
  expect_identical(one$alleles, two$alleles)
  expect_equal(one$records, two$records)
  expect_error(evolve(m, drift(1), n = 61), "even")
  expect_error(viability_selection(mates = 1, offspring_per_mate = 1), "at least 2")
})
