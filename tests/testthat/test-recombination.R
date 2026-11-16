# haplotype of locus l: paternal (column 2l - 1) or maternal (2l) allele
paternal <- function(g, l) g[, 2 * l - 1, ]
maternal <- function(g, l) g[, 2 * l, ]

test_that("genetic maps resolve to Haldane recombination fractions", {
  m <- evomod_model(additive(6, 1), mutation(0, 0.1), recombination = genetic_map(2, length = 0.5))
  r <- m$recombination
  expect_equal(r$type, "map")
  expect_equal(r$chromosome, c(1, 1, 1, 2, 2, 2))
  d <- 0.5 / 2
  expect_equal(r$rates, c(rep((1 - exp(-2 * d)) / 2, 2), 0.5, rep((1 - exp(-2 * d)) / 2, 2)))
  expect_equal(evomod_model(additive(4, 1), mutation(0, 0.1),
                            recombination = linkage(c(0.1, 0.5, 0)))$recombination$rates,
               c(0.1, 0.5, 0))
  expect_error(evomod_model(additive(4, 1), mutation(0, 0.1), recombination = linkage(c(0.1, 0.2))),
               "3")
  expect_error(linkage(0.7), "between 0 and 0.5")
  expect_error(evomod_model(additive(2, 1), mutation(0, 0.1), recombination = genetic_map(3)),
               "chromosomes")
})

# n individuals, one trait, all with paternal haplotype (1, 1, ...) and
# maternal haplotype (0, 0, ...): complete linkage disequilibrium
coupled <- function(loci, n) {
  g <- array(0, c(1, 2 * loci, n))
  g[, seq(1, 2 * loci, 2), ] <- 1
  g
}

test_that("crossovers happen with the given recombination fractions", {
  rates <- c(0.05, 0.3, 0.5, 0)
  m <- evomod_model(additive(5, 1), mutation(0, 0.1), recombination = linkage(rates))
  pop <- evolve(m, drift(1), n = 20000, initial = coupled(5, 20000), seed = 1)
  g <- genomes(pop)
  for (l in 1:4) {
    switched <- mean(c(paternal(g, l) != paternal(g, l + 1), maternal(g, l) != maternal(g, l + 1)))
    expect_lt(abs(switched - rates[l]), 0.01)
  }
})

test_that("linkage disequilibrium decays as (1 - r) per generation", {
  r <- 0.1
  n <- 5000
  m <- evomod_model(additive(2, 1), mutation(0, 0.1), recombination = linkage(r))
  pop <- evolve(m, drift(1), n = n, initial = coupled(2, n), seed = 2)
  D <- function(g) {
    a <- c(paternal(g, 1), maternal(g, 1)); b <- c(paternal(g, 2), maternal(g, 2))
    mean(a * b) - mean(a) * mean(b)
  }
  # every parent is a double heterozygote 11/00 (not random-mating genotype
  # frequencies), so the first generation gives D = 0.25 (1 - 2r); from then
  # on mating is random and D shrinks by (1 - r) per generation
  d1 <- D(genomes(pop))
  for (t in 2:10) pop <- evolve(pop, drift(1))
  expect_equal(d1, 0.25 * (1 - 2 * r), tolerance = 0.03)
  expect_equal(D(genomes(pop)), 0.25 * (1 - 2 * r) * (1 - r)^9, tolerance = 0.1)
})

test_that("without recombination whole haplotypes are inherited", {
  m <- evomod_model(additive(6, 2), mutation(0, 0.1), init_sd = 1,
                    recombination = linkage(rep(0, 5)))
  pop <- evolve(m, drift(1), n = 50, seed = 3)
  before <- genomes(pop)
  haps <- function(g) c(apply(g[, seq(1, 12, 2), ], 3, paste, collapse = ","),
                        apply(g[, seq(2, 12, 2), ], 3, paste, collapse = ","))
  pop <- evolve(pop, drift(3))
  expect_true(all(haps(genomes(pop)) %in% haps(before)))
})

test_that("free recombination is the default and a map with r = 0.5 is equivalent", {
  expect_equal(evomod_model(additive(3, 1), mutation(0, 0.1))$recombination$type, "free")
  m <- evomod_model(additive(5, 1), mutation(0, 0.1), recombination = linkage(rep(0.5, 4)))
  pop <- evolve(m, drift(1), n = 20000, initial = coupled(5, 20000), seed = 4)
  g <- genomes(pop)
  expect_lt(abs(mean(paternal(g, 1) != paternal(g, 5)) - 0.5), 0.01)
})

test_that("initial genomes are checked", {
  m <- evomod_model(additive(2, 1), mutation(0, 0.1))
  expect_error(evolve(m, drift(1), n = 10, initial = array(0, c(1, 4, 9))), "1 x 4 x 10")
})
