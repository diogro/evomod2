fl_model <- function(traits = 4) {
  evomod_model(additive(5, traits), mutation(0.01, 0.1), env_sd = 0.5, init_sd = 0.1)
}
theta_matrix <- function(pop) matrix(trait_means(pop)$theta, nrow = pop$model$traits)

test_that("the optimum jumps by step after every `every` generations", {
  pop <- evolve(fl_model(), fluctuating(200, every = 10, step = 0.2, omega = 1, bounds = NULL),
                n = 20, record_every = 1, seed = 1)
  th <- theta_matrix(pop)                     # traits x generations 0..200
  # changes[i] = G when the optimum differs between generations G - 1 and G;
  # jumps happen before generations 11, 21, ... (after 10, 20, ... generations)
  changes <- which(rowSums(abs(diff(t(th)))) > 0)
  expect_true(all(changes %% 10 == 1))
  d <- abs(diff(t(th)))
  expect_true(all(abs(d[d > 0] - 0.2) < 1e-12))
  expect_equal(length(changes), 19)           # unbounded: every jump is taken
})

test_that("correlated, independent and block structures", {
  corr <- theta_matrix(evolve(fl_model(), fluctuating(500, 10, 0.2, omega = 1, bounds = NULL),
                              n = 20, seed = 2))
  expect_true(all(apply(corr, 2, function(x) all(abs(x - x[1]) < 1e-12))))
  ind <- theta_matrix(evolve(fl_model(), fluctuating(5000, 10, 0.2, omega = 1, bounds = NULL,
                                                     structure = "independent"),
                             n = 20, record_every = 10, seed = 3))
  steps <- diff(t(ind))
  expect_lt(abs(cor(steps)[1, 2]), 0.15)
  blocks <- theta_matrix(evolve(fl_model(), fluctuating(5000, 10, 0.2, omega = 1, bounds = NULL,
                                                        structure = list(1:2, 3:4)),
                                n = 20, record_every = 10, seed = 4))
  expect_true(all(abs(blocks[1, ] - blocks[2, ]) < 1e-12))
  expect_true(all(abs(blocks[3, ] - blocks[4, ]) < 1e-12))
  s <- diff(t(blocks))
  expect_lt(abs(cor(s[, 1], s[, 3])), 0.15)
  # traits outside the blocks jump alone
  one <- theta_matrix(evolve(fl_model(3), fluctuating(3000, 10, 0.2, omega = 1, bounds = NULL,
                                                      structure = list(1:2)),
                             n = 20, record_every = 10, seed = 5))
  expect_true(all(abs(one[1, ] - one[2, ]) < 1e-12))
  expect_false(all(abs(one[1, ] - one[3, ]) < 1e-12))
})

test_that("bounds keep the optimum inside and make it stay half of the time at a bound", {
  pop <- evolve(fl_model(), fluctuating(20000, 10, 0.2, omega = 1, bounds = c(-1, 1),
                                        structure = "independent"),
                n = 10, record_every = 10, seed = 6)
  th <- theta_matrix(pop)
  expect_true(all(th >= -1 - 1e-9 & th <= 1 + 1e-9))
  # jump outcomes from the upper bound: stay or move down, 1/2 each
  x <- th[1, ]
  at_top <- which(abs(x[-length(x)] - 1) < 1e-9)
  stay <- mean(abs(x[at_top + 1] - 1) < 1e-9)
  expect_equal(stay, 0.5, tolerance = 0.15)
})

test_that("unbounded fluctuation is a random walk", {
  ends <- replicate(300, {
    pop <- evolve(fl_model(1), fluctuating(100, 10, 0.2, omega = 1, bounds = NULL),
                  n = 5, record_every = 100)
    pop$theta
  })
  # 9 jumps (after generations 10, 20, ..., 90) of +-0.2
  expect_equal(var(ends), 9 * 0.04, tolerance = 0.2)
  expect_equal(mean(ends), 0, tolerance = 0.1)
})

test_that("the optimum carries over between regimes and theta resets it", {
  m <- fl_model()
  pop <- evolve(m, list(fluctuating(100, 10, 0.2, omega = 1),
                        stabilizing(20, omega = 1)), n = 20, record_every = 1, seed = 7)
  th <- theta_matrix(pop)
  expect_equal(th[, 101], th[, 121])            # stabilizing keeps the last optimum
  expect_equal(pop$theta, th[, 121])
  pop <- evolve(pop, directional(10, delta_theta = rep(0.1, 4), omega = 1))
  expect_equal(pop$theta, th[, 121] + 1)
  pop <- evolve(pop, stabilizing(5, omega = 1, theta = rep(3, 4)))
  expect_equal(pop$theta, rep(3, 4))
  expect_error(evolve(pop, fluctuating(10, 5, 0.2, omega = 1)), "outside `bounds`")
})

test_that("fluctuating runs are reproducible and respond to selection", {
  m <- fl_model()
  reg <- fluctuating(300, 10, 0.2, omega = 1, structure = list(1:2, 3:4))
  a <- evolve(m, reg, n = 100, seed = 8)
  b <- evolve(m, reg, n = 100, seed = 8)
  expect_identical(a$records, b$records)
  expect_error(fluctuating(10, 5, 0.2, omega = 1, structure = "blocks"), "structure")
  expect_error(evolve(m, fluctuating(10, 5, 0.2, omega = 1, structure = list(1:2, 2:3)), n = 5),
               "disjoint")
  expect_error(fluctuating(10, 5, 0.2, omega = 1, bounds = c(1, -1)), "lower < upper")
})
