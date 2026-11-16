# Testing the claims of do O & Whitlock (2023) across genetic architectures.
#
# do O, I., & Whitlock, M. C. (2023). The evolution of genetic covariance and
# modularity as a result of multigenerational environmental fluctuation.
# Evolution Letters, 7(6), 457-466.
#
# Their claims, from simulations with a contribution matrix (z = C a):
#   1. correlated fluctuation of the optima gives more genetic covariance than
#      independent fluctuation;
#   2. most of the difference is pleiotropic, not linkage disequilibrium;
#   3. block-correlated fluctuation gives modularity;
#   4. the mutational covariance evolves to match the environment;
#   5. populations do best in the kind of environment they evolved in.
#
# The same environments are run on four architectures that differ in whether
# pleiotropy can evolve:
#   contribution_matrix  their model: continuous, evolvable pleiotropy
#   binary_b             evolvable pleiotropy through a binary B matrix
#   multilinear          fixed random epistasis: mutational effects depend on
#                        the genetic background, pleiotropic pattern fixed
#   additive             isotropic, fixed mutational covariance (control)
#
# Scales follow do O & Whitlock: 4 traits, N = 1000, 5000 generations, the
# optimum jumps by 0.2 every 10 generations within [-1, 1], omega = I, no
# environmental noise, 0.01 mutations per individual per generation of size
# 0.5 (the standard deviation used in their code). Statistics are averaged
# over the last 500 generations.
#
# Usage, one replicate per process:
#   Rscript fluctuating_environments.R simulate <architecture> <environment> <replicate> <outdir>
#   Rscript fluctuating_environments.R summarise <outdir> <summary.rds>
# architecture: contribution_matrix, binary_b, multilinear, additive
# environment:  correlated, independent, blocks

library(evomod2)

# ---- Settings ------------------------------------------------------------------

TRAITS <- 4
N <- 1000
GENERATIONS <- 5000
LAST <- 500                      # statistics over the last 500 generations
RECORD_EVERY <- 10
EVERY <- 10
STEP <- 0.2
BOUNDS <- c(-1, 1)
OMEGA <- diag(TRAITS)
BLOCKS <- list(1:2, 3:4)
SD <- 0.5                        # mutational effect size
U <- 0.01                        # mutations per individual per generation
LOCI <- 20                       # comparable to their 20 genes
TRANSPLANT <- 2000               # generations in each test environment
SHUFFLE_POINTS <- 10
ENVIRONMENTS <- c("correlated", "independent", "blocks")
ARCHITECTURES <- c("contribution_matrix", "binary_b", "multilinear", "additive")

structure_of <- function(env) switch(env, correlated = "correlated",
                                     independent = "independent", blocks = BLOCKS)

# models; `frozen` switches off the evolution of pleiotropy (C or B mutations)
make_model <- function(arch, frozen = FALSE, epsilon = NULL) {
  alleles <- 2 * LOCI
  switch(arch,
    contribution_matrix = evomod_model(
      contribution_matrix(TRAITS),                              # C = I, a = 1
      contribution_mutation(rate_weight = U / 2 / (2 * TRAITS),
                            rate_contribution = if (frozen) 0 else U / 2 / (2 * TRAITS^2),
                            sd = SD)),
    binary_b = evomod_model(
      binary_b(LOCI, TRAITS, init_b = 0.25),
      binary_b_mutation(rate = U / 2 / alleles, sd = SD,
                        rate_b = if (frozen) 0 else U / 2 / (alleles * TRAITS))),
    multilinear = evomod_model(multilinear(LOCI, TRAITS, epsilon),
                               mutation(rate = U / alleles, cov = SD)),
    additive = evomod_model(additive(LOCI, TRAITS), mutation(rate = U / alleles, cov = SD)))
}

# ---- Statistics ------------------------------------------------------------------

pairs_in <- function(m, which) {
  same <- outer(rep(seq_along(BLOCKS), lengths(BLOCKS)), rep(seq_along(BLOCKS), lengths(BLOCKS)), "==")
  off <- upper.tri(m)
  switch(which, all = m[off], within = m[off & same], between = m[off & !same])
}

covariance_stats <- function(m, prefix) {
  r <- cov2cor(m)
  out <- c(var = mean(diag(m)), cov = mean(pairs_in(m, "all")),
           cov_within = mean(pairs_in(m, "within")), cov_between = mean(pairs_in(m, "between")),
           cor = mean(pairs_in(r, "all")), cor_within = mean(pairs_in(r, "within")),
           cor_between = mean(pairs_in(r, "between")))
  setNames(out, paste0(prefix, "_", names(out)))
}

# ---- One replicate ------------------------------------------------------------------

simulate_replicate <- function(arch, env, replicate) {
  seed <- 1e6 * match(arch, ARCHITECTURES) + 1e5 * match(env, ENVIRONMENTS) + replicate
  set.seed(seed)
  epsilon <- if (arch == "multilinear") random_epsilon(TRAITS, sd = 0.1) else NULL
  model <- make_model(arch, epsilon = epsilon)
  regime <- function(g, bounds = BOUNDS)
    fluctuating(g, every = EVERY, step = STEP, omega = OMEGA, structure = structure_of(env),
                bounds = bounds)

  pop <- evolve(model, regime(GENERATIONS - LAST), n = N, record_every = RECORD_EVERY)
  # last generations in chunks, decomposing G after each
  decomp <- list()
  chunk <- LAST / SHUFFLE_POINTS
  for (i in seq_len(SHUFFLE_POINTS)) {
    pop <- evolve(pop, regime(chunk))
    decomp[[i]] <- g_decomposition(pop, method = "shuffle", shuffles = 5)
  }
  G <- g_matrices(pop)
  G_last <- Reduce(`+`, G[as.integer(names(G)) > GENERATIONS - LAST]) /
    sum(as.integer(names(G)) > GENERATIONS - LAST)
  pleio <- Reduce(`+`, lapply(decomp, `[[`, "genic")) / SHUFFLE_POINTS
  ld <- Reduce(`+`, lapply(decomp, `[[`, "ld")) / SHUFFLE_POINTS
  M <- mutational_covariance(pop, n = 20000)

  # transplant: pleiotropy frozen, unbounded optimum, each environment
  frozen <- make_model(arch, frozen = TRUE, epsilon = epsilon)
  fitness <- sapply(ENVIRONMENTS, function(test) {
    p <- pop
    p$model <- frozen
    p <- evolve(p, fluctuating(TRANSPLANT, every = EVERY, step = STEP, omega = OMEGA,
                               structure = structure_of(test), bounds = NULL,
                               theta = rep(0, TRAITS)))
    w <- mean_fitness(p)
    mean(w$mean_fitness[w$generation > max(w$generation) - TRANSPLANT / 10])
  })

  data.frame(architecture = arch, environment = env, replicate = replicate,
             t(covariance_stats(G_last, "G")), t(covariance_stats(pleio, "pleio")),
             t(covariance_stats(ld, "ld")[c("ld_var", "ld_cov", "ld_cov_within", "ld_cov_between")]),
             t(covariance_stats(M$per_mutation, "M")),
             t(setNames(fitness, paste0("fitness_", ENVIRONMENTS))))
}

save_replicate <- function(arch, env, replicate, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(simulate_replicate(arch, env, replicate),
          file.path(outdir, sprintf("%s_%s_%03d.rds", arch, env, replicate)))
}

summarise_runs <- function(outdir, file) {
  res <- do.call(rbind, lapply(list.files(outdir, "\\.rds$", full.names = TRUE), readRDS))
  saveRDS(list(results = res,
               settings = list(traits = TRAITS, n = N, generations = GENERATIONS, last = LAST,
                               every = EVERY, step = STEP, bounds = BOUNDS, sd = SD,
                               mutations_per_individual = U, loci = LOCI,
                               transplant = TRANSPLANT,
                               evomod2 = as.character(utils::packageVersion("evomod2")))),
          file, compress = "xz")
}

# ---- Command line -----------------------------------------------------------------

if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  cmd <- if (length(args)) args[1] else ""
  switch(cmd,
    simulate = save_replicate(args[2], args[3], as.integer(args[4]), args[5]),
    summarise = summarise_runs(args[2], args[3]),
    stop("usage: fluctuating_environments.R simulate <architecture> <environment> <replicate> <outdir> | summarise <outdir> <file>"))
}
