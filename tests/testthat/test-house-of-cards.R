test_that("house of cards mutations forget the old effect", {
  set.seed(1)
  # additive: alleles at 100, every allele mutates
  m <- evomod_model(additive(500, 2), hoc_mutation(rate = 1, cov = 0.5))
  out <- evomod2:::mutate_cpp(matrix(100, 2, 1000), m$mutation)
  changed <- out[, colSums(out != 100) > 0]
  expect_lt(abs(mean(changed)), 0.05)
  expect_equal(sd(as.vector(changed)), 0.5, tolerance = 0.05)
  # binary B: y is replaced, B untouched by y mutations
  bb <- binary_b_mutation(rate = 1, sd = 0.3, rate_b = 0, persistence = 0)
  g <- rbind(rep(100, 1000), matrix(1, 2, 1000))
  out <- evomod2:::mutate_cpp(g, bb)
  y <- out[1, out[1, ] != 100]
  expect_lt(abs(mean(y)), 0.05)
  expect_equal(sd(y), 0.3, tolerance = 0.05)
  expect_true(all(out[-1, ] == 1))
  # contribution matrix genes
  cm <- evomod_model(contribution_matrix(2), contribution_mutation(1, 1, 0.2, persistence = 0))
  out <- evomod2:::mutate_cpp(matrix(100, 1, 12), cm$mutation)
  expect_true(all(abs(out[out != 100]) < 2))
})

test_that("intermediate persistence regresses effects towards the reference", {
  set.seed(2)
  # low rate over many alleles, so that alleles are rarely hit twice
  m <- mutation(rate = 0.01, cov = 0.1, persistence = 0.5)
  mm <- evomod_model(additive(100000, 1), m)
  out <- evomod2:::mutate_cpp(matrix(10, 1, 200000), mm$mutation)
  changed <- out[out != 10]
  # median and MAD ignore the few alleles hit twice
  expect_equal(median(changed), 5, tolerance = 0.01)
  expect_equal(mad(changed), 0.1, tolerance = 0.1)
})

test_that("single house of cards mutations have effect: fresh draw minus old effect", {
  v <- 0.7; s <- 0.4
  m <- evomod_model(additive(10, 1), hoc_mutation(rate = 0.01, cov = s))
  pop <- evolve(m, drift(1), n = 20, seed = 3, initial = array(v, c(1, 20, 20)))
  pop$alleles[] <- v
  set.seed(4)
  mc <- mutational_covariance(pop, n = 40000)
  expect_equal(unname(mc$mean), -v, tolerance = 0.02)
  expect_equal(unname(mc$per_mutation[1, 1]), s^2, tolerance = 0.03)
})

test_that("neutral house of cards equilibrium: genic variance 2 L s^2 theta / (1 + theta)", {
  skip_on_cran()
  n <- 200; mu <- 0.00125; loci <- 100; s <- 0.5
  for (lc in list(wright_fisher(), moran())) {
    ne <- if (lc$type == "moran") n / 2 else n
    th <- 4 * ne * mu
    expected <- 2 * loci * s^2 * th / (1 + th)
    m <- evomod_model(additive(loci, 1), hoc_mutation(rate = mu, cov = s),
                      init_sd = s * sqrt(th / (1 + th)), life_cycle = lc)
    pop <- evolve(m, drift(3000), n = n, record_every = 100, seed = 5)
    g <- g_matrices(pop, "genic")
    v <- mean(sapply(g[as.integer(names(g)) > 500], function(x) x[1, 1]))
    expect_equal(v, expected, tolerance = 0.1)
  }
})

test_that("persistence is validated and shown", {
  expect_error(mutation(0.1, 0.1, persistence = 2), "between 0 and 1")
  expect_error(binary_b_mutation(0.1, 0.1, 0.1, persistence = -1), "between 0 and 1")
  expect_output(print(evomod_model(additive(2, 1), hoc_mutation(0.1, 0.1))), "house of cards")
  expect_output(print(evomod_model(additive(2, 1), mutation(0.1, 0.1))), "random walk")
  expect_equal(hoc_mutation(0.1, 0.2)$persistence, 0)
})
